import pandas as pd
from scipy.stats import mannwhitneyu

def filter_specificity(subnetwork_df, df_expression, target_ct, p_val_thresh=0.05, min_fc=1.2):
    """Validates target cell-type modular specificity against the expression dataset matrix."""
    if subnetwork_df.empty: return subnetwork_df
    
    mean_expr = df_expression.mean(axis=1)
    df_temp = pd.DataFrame({'expr': mean_expr})
    df_temp['Gene'] = df_temp.index.get_level_values("gene")
    df_temp['Cell'] = df_temp.index.get_level_values("celltype") + '_' + df_temp.index.get_level_values("patient")
    df_temp['CellType'] = df_temp.index.get_level_values("celltype")
    
    cell_expr_matrix = df_temp.pivot(index='Cell', columns='Gene', values='expr')
    cell_types = df_temp.drop_duplicates('Cell').set_index('Cell')['CellType']
    
    target_mask = (cell_types == target_ct)
    bg_mask = ~target_mask
    print(f"\nTesting specificity for '{target_ct}': {target_mask.sum()} target cells vs {bg_mask.sum()} background cells.")
    
    specific_subnetworks = []
    for _, row in subnetwork_df.iterrows():
        genes = [g for g in row['Genes'].split(', ') if g in cell_expr_matrix.columns]
        if not genes: continue
        
        cell_activity = cell_expr_matrix[genes].mean(axis=1)
        target_scores = cell_activity[target_mask]
        bg_scores = cell_activity[bg_mask]
        
        target_mean = target_scores.mean()
        bg_mean = bg_scores.mean()
        fold_change = (target_mean + 1e-9) / (bg_mean + 1e-9)
        
        if len(target_scores) < 3 or len(bg_scores) < 3:
            p_val = 0.0
        else:
            try:
                _, p_val = mannwhitneyu(target_scores, bg_scores, alternative='greater')
            except ValueError:
                p_val = 1.0
            
        if p_val <= p_val_thresh and fold_change >= min_fc:
            row_dict = row.to_dict()
            row_dict['Target_Mean'] = round(target_mean, 4)
            row_dict['Background_Mean'] = round(bg_mean, 4)
            row_dict['Fold_Change'] = round(fold_change, 3)
            specific_subnetworks.append(row_dict)
            
    final_df = pd.DataFrame(specific_subnetworks)
    if not final_df.empty:
        final_df = final_df.sort_values(by='Fold_Change', ascending=False)
    return final_df