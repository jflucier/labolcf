import csv
import os

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio import Entrez

def get_taxonomy_ids(accession, email="jean-francois.lucier@usherbrooke.ca"):
    """
    Retrieves taxonomy IDs for a list of NCBI accessions using BioPython's Entrez module.

    Args:
        accessions (list): A list of NCBI accession numbers (strings).
        email (str): Your email address (required by NCBI).

    Returns:
        list: A list of taxonomy IDs (strings), or None if an error occurs.
    """
    Entrez.email = email  # Always tell NCBI who you are!

    try:

        # Search for the accession in the Nucleotide database (or Protein, depending on your accessions)
        handle = Entrez.esearch(db="nucleotide", term=accession, retmode="xml") #Change to "protein" if needed.
        record = Entrez.read(handle)
        handle.close()

        if record["Count"] == '0':
            print(f"Accession {accession} not found.")
            return None

        gi_list = record["IdList"]

        if not gi_list:
            print(f"No GI number found for accession {accession}")
            return None

        gi = gi_list[0] #Take the first gi number.

        # Link the GI to the Taxonomy database
        link_handle = Entrez.elink(dbfrom="nucleotide", db="taxonomy", id=gi) #Change "nucleotide" to protein if needed.
        link_record = Entrez.read(link_handle)
        link_handle.close()

        if not link_record[0]["LinkSetDb"]:
            print(f"No taxonomy link found for accession {accession} (GI: {gi})")
            return None

        tax_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
        return tax_id

    except Exception as e:
        print(f"An error occurred: {e}")
        return None  # Return None if any errors occur.

# # Example usage:
# accession_list = ["CP000948.1", "NZ_CP011985.1", "NC_000913.3", "CP000107.1"] #Example accessions
# taxonomy_ids = get_taxonomy_ids(accession_list, email="your_email@example.com") #Replace with your email.

# if taxonomy_ids:
#     for accession, tax_id in zip(accession_list, taxonomy_ids):
#         print(f"Accession: {accession}, Taxonomy ID: {tax_id}")

filename = '/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/3_strains_sorted.tsv'  # Replace with your file name

data = []
row_count = 1
with open(filename, 'r', encoding='utf-8') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        # print(row)  # Process each row
        if row_count%10 == 0:
            print(f"# row parsed: {row_count}")

        row_count += 1

        if row[0] == '':
            continue

        if row[0] == '===':
            tax_id = get_taxonomy_ids(row[4])
            # genomes.append(row[4])
            # print(f"{row[4]}\t{tax_id}")
            data.append({'Accession': row[4], 'Taxonomy ID': tax_id if tax_id else None})  # Handle None return.

        # if row_count == 100:
        #     print(f"# row parsed: {row_count}")
        #     break

df = pd.DataFrame(data)
output_filename = "/storage/Documents/service/externe/lcfortier/20250226_abiF_genomic_context/3_strains_sorted.taxo.tsv"
df.to_csv(output_filename, sep='\t', index=False)
print(f"Results saved to {output_filename}")
# print("done")