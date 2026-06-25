import numpy as np
import pandas as pd
import fasttreeshap
from sklearn.tree import DecisionTreeClassifier
from joblib import Parallel, delayed

def _process_pattern(idx, X_dense_ct_slice, df_X, max_depth, target_ct):
    """Parallel kernel worker evaluating local trees and extracting positive explainer arrays."""
    y_counts = X_dense_ct_slice
    model = DecisionTreeClassifier(max_depth=max_depth, random_state=42)
    model.fit(df_X, y_counts)
    acc_score = model.score(df_X, y_counts)
    
    explainer = fasttreeshap.TreeExplainer(model)
    shap_values_raw = explainer.shap_values(df_X, check_additivity=False)
    
    if isinstance(shap_values_raw, list):
        shap_values_class_1 = shap_values_raw[1] 
    else:
        shap_values_class_1 = shap_values_raw[:, :, 1]
        
    mean_maxpos_shap = np.maximum(shap_values_class_1, 0).mean(axis=0)
    print(f"Pattern_{idx} - Total Count in {target_ct}: {y_counts.sum()} | Acc: {acc_score:.2f}")
    return idx, mean_maxpos_shap

def compute_shap_importance(X_binary, genes_vec, target_ct, max_depth=20, n_jobs=-1):
    """Fits localized tree networks in parallel to build the global gene-to-pattern mapping matrix."""
    print("Running Decision Trees and SHAP extraction in parallel...")
    df_X = pd.get_dummies(genes_vec, dtype=int)
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_pattern)(idx, X_binary[:, idx], df_X, max_depth, target_ct) 
        for idx in range(X_binary.shape[1])
    )
    
    shap_dict = {}
    for idx, shap_importance_array in results:
        shap_dict[idx] = shap_importance_array
        
    shap_df = pd.DataFrame(shap_dict, index=df_X.columns)
    print("SHAP extraction complete!")
    return shap_df.apply(pd.to_numeric, errors='coerce')