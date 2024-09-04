import argparse
import pandas as pd


def filter_blast(rep, out):
    pd_fs = pd.read_csv(rep, sep='\t', index_col="assembly")
    pd_fs_grouped = pd_fs.groupby("assembly")

    filtered_res = pd.DataFrame()
    for ass, rows in pd_fs_grouped:
        print(f"{ass}")
        rows.sort_values(by=['qcovs','pident'], ascending=False, inplace=True)
        best_qcov = 0
        best_pident = 0
        for index, row in rows.iterrows():
            curr_qcov = row["qcovs"]
            curr_pident = row["pident"]
            if curr_qcov >= best_qcov :
                if curr_pident >= best_pident:
                    tmp_df = pd.DataFrame([row])
                    filtered_res = pd.concat([filtered_res, tmp_df])
                    best_qcov = curr_qcov
                    best_pident = curr_pident
            else:
                break
    filtered_res.to_csv(
        out,
        sep="\t")

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-r",
        "--fs_report",
        help="foldseek output report with assembly label",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out",
        help="output filepath",
        required=True
    )

    args = argParser.parse_args()

    filter_blast(
        args.fs_report,
        args.out
    )
