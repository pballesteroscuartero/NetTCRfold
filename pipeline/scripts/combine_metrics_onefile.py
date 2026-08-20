import pandas as pd
import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folders containing all the models",
    )
    parser.add_argument(
        "-s",
        "--suffix",
    )
    return parser.parse_args()

def combine_metrics_onedb(path, db_name, file_name):
    contents = [f for f in path.iterdir() if f.is_dir()]
    df_db = pd.DataFrame()
    for f in contents:
        name = f.name
        df_folder = pd.read_csv(f'{f}/{file_name}')
        df_folder["db"] = db_name
        df_folder["param_comb"] = name
        df_db = pd.concat([df_db, df_folder], ignore_index=True)
    
    return df_db

def main():
    args=parse_args()
    path=Path(args.input)
    suffix = args.suffix
    df_smalldb_all = combine_metrics_onedb(path, "small_db", f"collected_af3metrics{suffix}.csv")
    df_smalldb_all.to_csv(path/ f"allresults_merged{suffix}.csv")

if __name__ == "__main__":
    main()