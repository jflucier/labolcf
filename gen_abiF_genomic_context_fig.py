import csv
import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def read_tsv_line_by_line(filename):
    """Reads a TSV file line by line and yields each row as a list."""
    try:
        with open(filename, 'r', encoding='utf-8') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                yield row
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def get_genepos_df(pre_pos, aft_pos):
    pre_pos.reverse()
    pre_data = []
    if len(pre_pos) > 5:
        print(f"#### length prepos > 5: {len(pre_pos)}")

    if len(aft_pos) > 5:
        print(f"#### length aft_pos > 5: {len(aft_pos)}")
        exit(0)

    for i in range(0,len(pre_pos)):
        pre_data.append([f"-{i+1}", pre_pos[i]])

    aft_data = []
    for i in range(0,len(aft_pos)):
        aft_data.append([f"{i+1}", aft_pos[i]])

    combined_data = pre_data + aft_data
    df = pd.DataFrame(combined_data, columns=["Pos", "Gene"])
    return df


# Example usage:
filename = '/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/3_strains_sorted.tsv'  # Replace with your file name
filter = '/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/abiF_coord.genomad_overlap.tsv'  # Replace with your file name

filter_df = pd.read_csv(filter, sep='\t', header=0)
valid_genomes = filter_df['abi_genome_acc']
# Make valid_genomes unique
unique_valid_genomes = valid_genomes.unique()

pre_pos = []
aft_pos = []
reached_abiF = False
all_dfs = []
row_count = 0
prev_genome_acc = ""
with open(filename, 'r', encoding='utf-8') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        # print(row)  # Process each row
        if row_count%1000 == 0:
            print(f"# row parsed: {row_count}")

        row_count += 1

        if row[0] == '':
            continue

        # if len(aft_pos) > 5:
        #     print(f"#### length aft_pos > 5: {len(aft_pos)}")
        #     exit(0)

        if row[0] == '===':
            if row_count == 1:
                parts = row[1].split("_")
                if len(parts[0]) == 2:
                    prev_genome_acc = parts[0] + "_" + parts[1]
                else:
                    prev_genome_acc = parts[0]
                continue

            if len(pre_pos) == 0 and len(aft_pos) == 0:
                reached_abiF = False
                pre_pos = []
                aft_pos = []
                continue

            if reached_abiF == False:
                print(f"Whats going on?: {row}")
                # skip for now
                reached_abiF = False
                pre_pos = []
                aft_pos = []
                continue

            if prev_genome_acc in unique_valid_genomes:
                partial_df = get_genepos_df(pre_pos, aft_pos)
                all_dfs.append(partial_df)

            # reset all var
            reached_abiF = False
            pre_pos = []
            aft_pos = []
            parts = row[1].split("_")
            if len(parts[0]) == 2:
                prev_genome_acc = parts[0] + "_" + parts[1]
            else:
                prev_genome_acc = parts[0]
            continue

        if (row[9] == 'Abi_2' or row[9] == 'AbiF') and reached_abiF == True :
            # weird case where abiF and abi 2 detected, reset aft_pos array to ignore in between
            aft_pos = []
            continue

        if row[9] == 'Abi_2' or row[9] == 'AbiF':
            reached_abiF = True
            continue

        if reached_abiF == False:
            pre_pos.append(row[8])
        else:
            aft_pos.append(row[8])

# add last item
partial_df = get_genepos_df(pre_pos, aft_pos)
all_dfs.append(partial_df)

if all_dfs: #check to see if the list is not empty.
    pos_df = pd.concat(all_dfs, ignore_index=True)
else:
    pos_df = pd.DataFrame(columns=["Pos", "Gene"])

gene_counts = pos_df.groupby(['Pos','Gene'])["Gene"].count()
gene_counts_df = gene_counts.reset_index(name='Count')
total_counts_per_pos = gene_counts_df.groupby('Pos')['Count'].transform('sum')
gene_counts_df['Proportion'] = gene_counts_df['Count'] / total_counts_per_pos

gene_counts_df.to_csv(
    f"/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/abiF_fig.category.abif_genomad.tsv",
    index=False,
    sep='\t'
)

print("done")