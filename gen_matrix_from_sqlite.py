import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sqlite3
import numpy as np

def gen_heat(mat, out):
    # data_pivoted = mat.pivot(index=["slpa_type","assembly"], columns=["pident"], values=["pident"])
    # ax = sns.heatmap(mat)
    f = plt.figure()
    # ax = sns.clustermap(mat, z_score=1, method="complete", metric="euclidean", cmap='RdBu', annot=True,
    #            annot_kws={"size": 7}, vmin=-2, vmax=2, figsize=(15,12))
    ax = sns.clustermap(mat, method="complete", metric="euclidean", cmap='RdBu', annot=True,
                        annot_kws={"size": 7}, vmin=-0, vmax=50, figsize=(15, 12))
    # plt.show()
    plt.tight_layout()
    plt.savefig(out)

def gen_matrix(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    rbps = pd.read_sql_query("select distinct qseqid from RBP_list_AA_2024_blastout_blastp_qcov95_ident90", con)
    slpas = pd.read_sql_query("select distinct slpa from Slpa_blast_qcov90_ident60", con)

    matrix = pd.DataFrame(0, index=rbps['qseqid'], columns=slpas['slpa'])

    for rbp in rbps['qseqid']:
        for slpa in slpas['slpa']:
            c = pd.read_sql_query(f"""
            select
              count(distinct b.assembly) c
            from RBP_list_AA_2024_blastout_blastp_qcov95_ident90 r
            JOIN Slpa_blast_qcov90_ident60 b on r.assembly=b.assembly
            where r.qseqid = '{rbp}' and b.slpa = '{slpa}'
            """, con)
            # print(f"rbp={rbp} slpa={slpa} val={c.iloc[0]['c']}")
            matrix.at[rbp, slpa] = c.iloc[0]['c']
    con.close()

    return matrix

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-db",
        help="sqlite file",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out",
        help="output filepath",
        required=True
    )

    args = argParser.parse_args()

    m = gen_matrix(
        args.db,
        args.out
    )

    gen_heat(m, args.out)