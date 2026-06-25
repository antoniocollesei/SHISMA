# shisma/borf.py

import os
import sys
import importlib.util
import numpy as np
import pandas as pd
from fast_borf import BorfBuilder
from fast_borf.pipeline.zero_columns_remover import ZeroColumnsRemover
from fast_borf.pipeline.reshaper import ReshapeTo2D
from fast_borf.pipeline.to_scipy import ToScipySparse

def _load_borf_config(config_path=None):
    """Loads BORF_CONFIG. Defaults to the borf_configuration.py sitting alongside this file."""
    if config_path is None:
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(pkg_dir, "borf_configuration.py")
        
    if not os.path.exists(config_path):
        print(f"Warning: Configuration file '{config_path}' not found.")
        print("Falling back to default backup window configuration matrix...")
        return [
            {"window_size": 2, "stride": 1, "dilation": 1, "word_length": 2, "alphabet_size": 3},
            {"window_size": 3, "stride": 1, "dilation": 1, "word_length": 3, "alphabet_size": 3},
            {"window_size": 4, "stride": 1, "dilation": 1, "word_length": 4, "alphabet_size": 3}
        ]
        
    spec = importlib.util.spec_from_file_location("dynamic_borf_config", config_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["dynamic_borf_config"] = module
    spec.loader.exec_module(module)
    return getattr(module, "BORF_CONFIG")

def run_borf_transformation(df_expression, ppi_df, target_ct, config_file_path=None, n_jobs=-1):
    """Processes array matrices and pipes variables through fast-borf features extraction maps."""
    print("1. Filtering for PPI genes & Structuring Multi-Index Slices...")
    ppi_genes = np.unique(ppi_df[['gene1', 'gene2']].values.astype(str))
    
    if not isinstance(df_expression.index, pd.MultiIndex):
        index_dims = df_expression.index.str.rsplit('_', n=2, expand=True)
        index_dims.names = ["gene", "patient", "celltype"] 
        df_expression.index = index_dims
    
    df_filtered = df_expression[df_expression.index.get_level_values("gene").isin(ppi_genes)]
    print(f" -> Reduced to {len(df_filtered.index.get_level_values('gene').unique())} valid PPI genes.")
    
    df_target_ct = df_filtered.xs(target_ct, level="celltype")
    genes_vec = df_target_ct.index.get_level_values("gene").values
    X_raw = df_target_ct.values[:, np.newaxis, :]
    
    borf_configurations_profile = _load_borf_config(config_file_path)
    
    print(f"\n2. Executing BORF Transformation pipeline...")
    builder = BorfBuilder(
        n_jobs=n_jobs, 
        configs=borf_configurations_profile,
        pipeline_objects=[
            (ReshapeTo2D, dict(keep_unraveled_index=True)), 
            (ZeroColumnsRemover, dict(axis=0)), 
            (ToScipySparse, dict())
        ]
    )
    
    borf_model = builder.build(X_raw)
    X_dense = borf_model.fit_transform(X_raw).toarray()
    X_binary = (X_dense != 0).astype(int)
    print(f"BORF Matrix Generation Complete. Shape: {X_binary.shape}")
    
    return X_binary, genes_vec, borf_model, X_raw, df_target_ct