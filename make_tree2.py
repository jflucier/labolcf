import sys
import argparse
from ete3 import NCBITaxa, TreeStyle

def plot_family_tree_from_ids(taxon_ids, output_image_file="family_tree.png"):
    """
    Builds and plots a phylogenetic tree where families related to the
    input list of NCBI taxon IDs are the leaves.

    Args:
        taxon_ids (list): A list of NCBI taxon IDs (integers). The families
                          related to these IDs will become the leaves.
        output_image_file (str): Path to save the tree image.
    """
    ncbi = NCBITaxa()

    try:
        valid_taxon_ids = []
        invalid_taxon_ids = []
        for tid_str in taxon_ids:
            tid_str = tid_str.strip()  # Remove leading/trailing whitespace
            if tid_str:  # Ignore empty lines
                try:
                    taxon_id = int(tid_str)
                    valid_taxon_ids.append(taxon_id)
                except ValueError:
                    invalid_taxon_ids.append(tid_str)

        if invalid_taxon_ids:
            print(f"Warning: The following are not valid integer taxon IDs and will be ignored: {invalid_taxon_ids}")

        if not valid_taxon_ids:
            print("Error: No valid integer taxon IDs provided.")
            return

        all_descendant_families = set()
        for taxon_id in valid_taxon_ids:
            try:
                descendant_families = ncbi.get_descendant_taxa(taxon_id, rank='family')
                all_descendant_families.update(descendant_families)
            except ValueError as e:
                print(f"Warning: Could not find taxon ID {taxon_id} in the NCBI database: {e}")

        if not all_descendant_families:
            print(f"No families found for the valid taxon IDs provided: {valid_taxon_ids}")
            return

        # Get the phylogenetic topology for these families
        tree = ncbi.get_topology(list(all_descendant_families))

        # Prune the tree to keep only the family nodes as leaves
        families_names = set(ncbi.get_taxid_translator(all_descendant_families).values())

        for node in tree.traverse("postorder"):
            if node.is_leaf() and node.name not in families_names:
                node.delete()
            elif not node.is_leaf():
                has_family_leaf = False
                for leaf in node.iter_leaves():
                    if leaf.name in families_names:
                        has_family_leaf = True
                        break
                if not has_family_leaf and node != tree:  # Don't remove the root
                    node.delete()

        # Create a TreeStyle
        def my_layout(node):
            node.img_style['size'] = 15

        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = False
        ts.branch_vertical_margin = 10
        ts.layout_fn = my_layout

        # Render and save the tree
        tree.render(output_image_file, tree_style=ts)
        print(f"Family-level tree saved to: {output_image_file}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Build and plot a phylogenetic tree with families as leaves based on a comma-separated list of NCBI taxon IDs.")
    # parser.add_argument("--ids", dest="taxon_ids_str", required=True, help="A comma-separated list of NCBI taxon IDs (integers).")
    # parser.add_argument("--out", dest="output_image", default="family_tree_from_list.png", help="Path to save the tree image (default: family_tree_from_list.png).")
    #
    # args = parser.parse_args()
    #
    # taxon_ids_list = args.taxon_ids_str.split(',')
    # plot_family_tree_from_ids(taxon_ids_list, args.output_image)

    # tax = "470,546,571,575,584,587,588,615,630,672,735,750,818,853,907,1245,1246,1247,1249,1254,1255,1260,1281,1283,1290,1292,1294,1296,1303,1304,1305,1307,1311,1313"
    # taxon_ids_list = tax.split(',')
    # out = "/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/test.png"
    # plot_family_tree_from_ids(taxon_ids_list, out)


