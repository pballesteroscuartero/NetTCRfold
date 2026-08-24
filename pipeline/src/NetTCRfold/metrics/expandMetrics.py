import pandas as pd
import numpy as np
import ast
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Expand metrics from chain pair matrices into separate columns")
    parser.add_argument("--input_csv", "-i", type=str, required=True, help="Path to the input CSV file containing the metrics")
    parser.add_argument("--output_csv", "-o", type=str, required=True, help="Path to save the expanded CSV file")
    return parser.parse_args()

def strictly_upper_triangle(mat):
    """Return the flattened strictly upper triangle of a square 2D list or array"""
    arr = np.array(mat)
    return arr[np.triu_indices_from(arr, k=1)].tolist()

def off_diagonal(mat):
    """Return all off-diagonal elements flattened (row-wise)"""
    arr = np.array(mat)
    return arr[~np.eye(arr.shape[0], dtype=bool)]

def flatten_dict(d, key_map, mode="upper"):
    if isinstance(d, str):
        try:
            d = ast.literal_eval(d)
        except Exception:
            d = eval(d.replace('nan', 'np.nan'))
    keys_in_row = [k for k in key_map.keys() if k in d]
    chain_order = [key_map[k] for k in keys_in_row]
    values = []
    key_pairs = []
    for i, r_key in enumerate(keys_in_row):
        for j, c_key in enumerate(keys_in_row):
            if mode =="upper":
                if j > i:
                    values.append(d[r_key][c_key])
                    key_pairs.append((key_map[r_key], key_map[c_key]))
            elif mode =="offdiag":
                if i != j:
                    values.append(d[r_key][c_key])
                    key_pairs.append((key_map[r_key], key_map[c_key]))

    return values, key_pairs

def extract_dictionary_into_df(df, row, idx, metrics, dict_key_map, mode_flatten = "offdiag"):
    #record = {}
    for metric in metrics:
        matrix = row[metric]
        vals, order =flatten_dict(matrix, dict_key_map, mode=mode_flatten)
        for col_name, val in zip(order, vals):
            colname = f"{col_name[0]}_{col_name[1]}_{metric}"
            df.at[idx, colname] = val

def expand_chain_pair_metrics_gen(df):
    #I think AF provides the metrics in alphabetical order of the chain, not in the chain provided order in json file
    dict_key_map = {'A': 'MHCA', 'B': 'MHCB', 'C': 'pep', 'D': 'TRA', 'E': 'TRB'}
    nchains = df["chain_iptm"].apply(ast.literal_eval).apply(lambda x: len(x))
    chain_names_list=[]
    for i in nchains:
        if i == 4:
            chain_names_list.append(['MHCA', 'pep', 'TRA', 'TRB'])
        elif i == 5:
            chain_names_list.append(['MHCA','MHCB', 'pep', 'TRA', 'TRB'])

    df["chain_names"] = chain_names_list

    for idx, row in df.iterrows():
        chains = row["chain_names"]
        for col_name, value_iptm, value_ptm in zip(row['chain_names'], ast.literal_eval(row['chain_iptm']), ast.literal_eval(row['chain_ptm'])):
            df.at[idx, f"{col_name}_iptm"] = value_iptm
            df.at[idx, f"{col_name}_ptm"] = value_ptm

        matrix_iptm = ast.literal_eval(row["chain_pair_iptm"])  
        upper_vals_iptm = strictly_upper_triangle(matrix_iptm) #Symmetric, so upper is good
        pair_iptm_columns_upper = [f'{c}_{r}_iptm' for i, r in enumerate(row["chain_names"]) for j, c in enumerate(row["chain_names"]) if j > i]
        for col_name, val in zip(pair_iptm_columns_upper, upper_vals_iptm):
            df.at[idx, col_name] = val
            
        matrix_pae = ast.literal_eval(row["chain_pair_pae_min"])
        vals_pae = off_diagonal(matrix_pae)
        pair_pae_columns = [f"{c}_{r}_pae_min" 
                            for i, r in enumerate(row["chain_names"])
                            for j, c in enumerate(row["chain_names"])
                            if i != j
                            ]
        for col_name, val in zip(pair_pae_columns, vals_pae):
            df.at[idx, col_name] = val

        metrics = ["cdr_metric_mean_chain",
                   "ipsae", "ipsae_d0chn", "ipsae_d0dom"]

        extract_dictionary_into_df(df, row, idx, metrics, dict_key_map, mode_flatten="offdiag")
    
    return df   


def main():
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    expanded_df = (
        df.groupby("param_comb", group_keys=False)
        .apply(expand_chain_pair_metrics_gen)
    )
    expanded_df = expanded_df.drop(columns = ["chain_iptm", "chain_ptm", "chain_pair_iptm",
                                               "chain_pair_pae_min", 
                                               "cdr_metric_mean_chain",
                                               "ipsae", "ipsae_d0chn", "ipsae_d0dom",
                                                "chain_names"])
    expanded_df = expanded_df.drop(columns=[col for col in expanded_df.columns if "MHCB" in col])
    expanded_df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()

