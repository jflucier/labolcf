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
    # ax = sns.clustermap(mat, method="complete", metric="euclidean", cmap='RdBu', annot=True,
    #                     annot_kws={"size": 7}, vmin=-0, vmax=25, figsize=(15, 12))
    # ax = sns.clustermap(mat, method="complete", metric="euclidean", cmap='RdBu', figsize=(15, 12))
    ax = sns.clustermap(mat, method="complete", metric="euclidean", cmap=['white','black'], figsize=(15, 12))
    # plt.show()
    plt.tight_layout()
    plt.savefig(out)

def gen_matrix_slpa_rbp(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    rbps = pd.read_sql_query("select distinct qseqid from RBP_list_AA_2024_blastout_blastp_qcov95_besthit", con)
    slpas = pd.read_sql_query("select distinct slpa from Slpa_blast_qcov90_ident60", con)

    matrix = pd.DataFrame(0, index=rbps['qseqid'], columns=slpas['slpa'])

    for rbp in rbps['qseqid']:
        for slpa in slpas['slpa']:
            c = pd.read_sql_query(f"""
            select
              count(distinct b.assembly) c
            from RBP_list_AA_2024_blastout_blastp_qcov95_besthit r
            JOIN Slpa_blast_qcov90_ident60 b on r.assembly=b.assembly
            where r.qseqid = '{rbp}' and b.slpa = '{slpa}'
            """, con)
            # print(f"rbp={rbp} slpa={slpa} val={c.iloc[0]['c']}")
            matrix.at[rbp, slpa] = c.iloc[0]['c']
    con.close()

    return matrix

def gen_matrix_assembly_slpa_rbp(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    assembly = pd.read_sql_query("select distinct assembly from Slpa_blast_qcov90_ident60", con)
    rbps = pd.read_sql_query("select distinct qseqid x from RBP_list_AA_2024_blastout_blastp_qcov95_besthit", con)
    slpas = pd.read_sql_query("select distinct slpa x from Slpa_blast_qcov90_ident60", con)
    all_cols = pd.concat([rbps,slpas],ignore_index=True)
    # all_cols.loc[len(all_cols)] = {'x': 'prophage'}
    # all_cols.loc[len(all_cols)] = {'x': 'diffocine'}

    matrix = pd.DataFrame(int(0), index=assembly['assembly'], columns=all_cols['x'])
    for col in matrix.columns:
        matrix[col].values[:] = int(0)

    # for ass in assembly['assembly']:
    for slpa in slpas['x']:
        c = pd.read_sql_query(f"""
        select
          b.assembly ass
        from Slpa_blast_qcov90_ident60 b
        where b.slpa = '{slpa}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], slpa] = int(1)

    for rbp in rbps['x']:
        c = pd.read_sql_query(f"""
        select
          r.assembly ass
        from Slpa_blast_qcov90_ident60 b 
        join RBP_list_AA_2024_blastout_blastp_qcov95_besthit r on b.assembly=r.assembly
        where r.qseqid = '{rbp}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], rbp] = int(1)

    # c = pd.read_sql_query(f"""
    # select
    #   r.ERR017367assembly ass
    # from Slpa_blast_qcov90_ident60 b
    # join genomad_prophage r on b.assembly=r.ERR017367assembly
    # """, con)
    #
    # # assign 1 to all assembly found in matrix
    # for index, row in c.iterrows():
    #     matrix.at[row["ass"], 'prophage'] = int(1)
    #
    # c = pd.read_sql_query(f"""
    # select
    #   r.ERR017367assembly ass
    # from Slpa_blast_qcov90_ident60 b
    # join genomad_diffocine r on b.assembly=r.ERR017367assembly
    # """, con)
    #
    # # assign 1 to all assembly found in matrix
    # for index, row in c.iterrows():
    #     matrix.at[row["ass"], 'diffocine'] = int(1)

    con.close()

    matrix.to_csv("/storage/Documents/service/externe/lcfortier/20240618_cdiff_assembly_analysis/20240618_cdiff_assembly_analysis.slpa_vs_rbp.tsv")

    return matrix

def gen_matrix_assembly_slpa_rbp_genomad(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    assembly = pd.read_sql_query("select distinct assembly from Slpa_blast_qcov90_ident60", con)
    rbps = pd.read_sql_query("select distinct qseqid x from RBP_list_AA_2024_blastout_blastp_qcov95_besthit", con)
    slpas = pd.read_sql_query("select distinct slpa x from Slpa_blast_qcov90_ident60", con)
    all_cols = pd.concat([rbps,slpas],ignore_index=True)
    all_cols.loc[len(all_cols)] = {'x': 'prophage'}
    all_cols.loc[len(all_cols)] = {'x': 'diffocine'}

    matrix = pd.DataFrame(int(0), index=assembly['assembly'], columns=all_cols['x'])
    for col in matrix.columns:
        matrix[col].values[:] = int(0)

    # for ass in assembly['assembly']:
    for slpa in slpas['x']:
        c = pd.read_sql_query(f"""
        select
          b.assembly ass
        from Slpa_blast_qcov90_ident60 b
        where b.slpa = '{slpa}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], slpa] = int(1)

    for rbp in rbps['x']:
        c = pd.read_sql_query(f"""
        select
          r.assembly ass
        from Slpa_blast_qcov90_ident60 b 
        join RBP_list_AA_2024_blastout_blastp_qcov95_besthit r on b.assembly=r.assembly
        where r.qseqid = '{rbp}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], rbp] = int(1)

    c = pd.read_sql_query(f"""
    select
      r.ERR017367assembly ass
    from Slpa_blast_qcov90_ident60 b
    join genomad_prophage r on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'prophage'] = int(1)

    c = pd.read_sql_query(f"""
    select
      r.ERR017367assembly ass
    from Slpa_blast_qcov90_ident60 b
    join genomad_diffocine r on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'diffocine'] = int(1)

    con.close()

    matrix.to_csv("/storage/Documents/service/externe/lcfortier/20240618_cdiff_assembly_analysis/20240618_cdiff_assembly_analysis.slpa_vs_rbp.tsv")

    return matrix

def gen_matrix_assembly_slpa_rbp_genomad_prophage_filter(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    assembly = pd.read_sql_query('''select 
      r.ERR017367assembly assembly
    from genomad_prophage r
    join Slpa_blast_qcov90_ident60 b on b.assembly=r.ERR017367assembly''', con)
    rbps = pd.read_sql_query("select distinct qseqid x from RBP_list_AA_2024_blastout_blastp_qcov95_besthit", con)
    slpas = pd.read_sql_query("select distinct slpa x from Slpa_blast_qcov90_ident60", con)
    all_cols = pd.concat([rbps,slpas],ignore_index=True)
    all_cols.loc[len(all_cols)] = {'x': 'prophage'}
    all_cols.loc[len(all_cols)] = {'x': 'diffocine'}

    matrix = pd.DataFrame(int(0), index=assembly['assembly'], columns=all_cols['x'])
    for col in matrix.columns:
        matrix[col].values[:] = int(0)

    # for ass in assembly['assembly']:
    for slpa in slpas['x']:
        c = pd.read_sql_query(f"""
        select
          r.ERR017367assembly ass
        from genomad_prophage r
        join Slpa_blast_qcov90_ident60 b on b.assembly=r.ERR017367assembly
        where b.slpa = '{slpa}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], slpa] = int(1)

    for rbp in rbps['x']:
        c = pd.read_sql_query(f"""
        select
          p.ERR017367assembly ass
        from genomad_prophage p
        join RBP_list_AA_2024_blastout_blastp_qcov95_besthit r on b.assembly=r.assembly
        join Slpa_blast_qcov90_ident60 b on b.assembly=p.ERR017367assembly
        where r.qseqid = '{rbp}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], rbp] = int(1)

    c = pd.read_sql_query(f"""
    select
      r.ERR017367assembly ass
    from genomad_prophage r
    join Slpa_blast_qcov90_ident60 b on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'prophage'] = int(1)

    c = pd.read_sql_query(f"""
    select
      p.ERR017367assembly ass
    from genomad_prophage r
    join genomad_diffocine p on p.ERR017367assembly=r.ERR017367assembly    
    join Slpa_blast_qcov90_ident60 b on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'diffocine'] = int(1)

    con.close()

    matrix.to_csv("/storage/Documents/service/externe/lcfortier/20240618_cdiff_assembly_analysis/20240618_cdiff_assembly_analysis.slpa_vs_rbp.tsv")

    return matrix

def gen_matrix_assembly_rbp_defensefinder(sqlitedb, out):
    con = sqlite3.connect(sqlitedb)

    # Load the data into a DataFrame
    assembly = pd.read_sql_query("select distinct assembly from Slpa_blast_qcov90_ident60", con)
    rbps = pd.read_sql_query("select distinct r.qseqid x from RBP_list_AA_2024_blastout_blastp_qcov95_besthit r join Slpa_blast_qcov90_ident60 b on b.assembly=r.assembly", con)
    deffinder = pd.read_sql_query("select distinct r.type x from defensefinder_results r join Slpa_blast_qcov90_ident60 b on b.assembly=r.assembly", con)
    all_cols = pd.concat([rbps,deffinder],ignore_index=True)
    all_cols.loc[len(all_cols)] = {'x': 'prophage'}
    all_cols.loc[len(all_cols)] = {'x': 'diffocine'}

    matrix = pd.DataFrame(int(0), index=assembly['assembly'], columns=all_cols['x'])
    for col in matrix.columns:
        matrix[col].values[:] = int(0)

    # for ass in assembly['assembly']:
    for type in deffinder['x']:
        c = pd.read_sql_query(f"""
        select
          b.assembly ass
        from Slpa_blast_qcov90_ident60 b 
        join defensefinder_results r on b.assembly=r.assembly
        where r.type = '{type}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], type] = int(1)

    for rbp in rbps['x']:
        c = pd.read_sql_query(f"""
        select
          r.assembly ass
        from Slpa_blast_qcov90_ident60 b 
        join RBP_list_AA_2024_blastout_blastp_qcov95_besthit r on b.assembly=r.assembly
        where r.qseqid = '{rbp}'
        """, con)

        # assign 1 to all assembly found in matrix
        for index, row in c.iterrows():
            matrix.at[row["ass"], rbp] = int(1)

    c = pd.read_sql_query(f"""
    select
      r.ERR017367assembly ass
    from Slpa_blast_qcov90_ident60 b
    join genomad_prophage r on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'prophage'] = int(1)

    c = pd.read_sql_query(f"""
    select
      r.ERR017367assembly ass
    from Slpa_blast_qcov90_ident60 b
    join genomad_diffocine r on b.assembly=r.ERR017367assembly
    """, con)

    # assign 1 to all assembly found in matrix
    for index, row in c.iterrows():
        matrix.at[row["ass"], 'diffocine'] = int(1)

    con.close()

    matrix.to_csv("/storage/Documents/service/externe/lcfortier/20240618_cdiff_assembly_analysis/20240618_cdiff_assembly_analysis.slpa_vs_rbp.tsv")

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

    m = gen_matrix_assembly_rbp_defensefinder(
        args.db,
        args.out
    )

    gen_heat(m, args.out)