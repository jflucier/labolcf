from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import pandas as pd

def create_tree_from_mash_dist(mash_dist_file, output_newick_file):
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

        print(f"Phylogenetic tree saved to: {output_newick_file}")

    except FileNotFoundError:
        print(f"Error: File not found: {mash_dist_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:
mash_dist_file = "/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/abiF_genomad_overlap.acc.distances.txt"  # Replace with your Mash dist output file
output_newick_file = "/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/abiF_genomad_overlap.acc.distances.nwk" # Replace with your desired output path

create_tree_from_mash_dist(mash_dist_file, output_newick_file)