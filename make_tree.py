import argparse

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import pandas as pd
from ete3 import Tree, TreeStyle
from matplotlib import pyplot as plt


def create_tree_from_mash_dist(mash_dist_file, output_newick_file,output_png_file):
    """
    Creates a phylogenetic tree from Mash dist output.

    Args:
        mash_dist_file (str): Path to the Mash dist output file.
        output_newick_file (str): Path to save the Newick tree file.
    """

    try:
        # Read the mash dist file into a pandas dataframe
        df = pd.read_csv(mash_dist_file, sep='\t', header=None, names=['Genome1', 'Genome2', 'Distance', 'p-value', 'shared-hashes'])

        # Extract unique genome names
        genome_names = sorted(list(set(df['Genome1'].unique()) | set(df['Genome2'].unique())))

        # Create a lower triangular distance matrix
        distance_matrix = [[0.0] * i for i in range(1, len(genome_names))]

        # Populate the lower triangular distance matrix
        name_to_index = {name: i for i, name in enumerate(genome_names)}
        for _, row in df.iterrows():
            i = name_to_index[row['Genome1']]
            j = name_to_index[row['Genome2']]
            if i > j:
                distance_matrix[i - 1][j] = row['Distance']
            elif j > i:
                distance_matrix[j - 1][i] = row['Distance']

        # Adjust names to match the lower triangle matrix
        adjusted_names = genome_names[1:]

        # Create a DistanceMatrix object
        dist_matrix = DistanceMatrix(adjusted_names, distance_matrix)

        # Build the Neighbor-Joining tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dist_matrix)

        # Save the tree to a Newick file
        Phylo.write(tree, output_newick_file, "newick")

        fig = plt.figure(figsize=(10, 20), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        plt.show()
        plt.savefig(output_png_file)

        print(f"Phylogenetic tree saved to: {output_newick_file}")

    except FileNotFoundError:
        print(f"Error: File not found: {mash_dist_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

# def plot_newick_ete3(newick_file, output_image_file):
#     """
#     Plots a Newick tree using ete3.
#
#     Args:
#         newick_file (str): Path to the Newick tree file.
#         output_image_file (str): Path to save the image (e.g., PNG, SVG).
#     """
#
#     try:
#         t = Tree(newick_file)
#         ts = TreeStyle()
#         ts.show_leaf_name = True  # Show leaf names
#         ts.show_branch_length = True # show branch length
#
#         # Customize the tree style (optional)
#         ts.branch_vertical_margin = 10
#         ts.scale = 100
#
#         t.render(output_image_file, tree_style=ts)
#         print(f"Tree plot saved to: {output_image_file}")
#
#     except FileNotFoundError:
#         print(f"Error: File not found: {newick_file}")
#     except Exception as e:
#         print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a phylogenetic tree from Mash dist output.")
    parser.add_argument("--in", dest="mash_dist_file", required=True, help="Path to the Mash dist output file.")
    parser.add_argument("--out", dest="output_newick_file", required=True, help="Path to save the Newick tree file.")
    parser.add_argument("--png", dest="output_png_file", required=True, help="Path to save the tree PNG file.")

    args = parser.parse_args()

    create_tree_from_mash_dist(
        args.mash_dist_file,
        args.output_newick_file,
        args.output_png_file
    )

