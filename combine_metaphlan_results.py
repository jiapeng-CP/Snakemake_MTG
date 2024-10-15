import os
import pandas as pd
from glob import glob

def main():
    # Path to the directory containing MetaPhlAn result files
    input_directory = "./MetaPhlan/*.txt"  # Adjust this path as necessary
    output_file = "combined_metaphlan_results.tsv"

    # Initialize an empty list to store data frames
    dfs = []

    # List all MetaPhlAn output files
    metaphlan_files = glob(input_directory)

    if not metaphlan_files:
        print("No MetaPhlAn result files found in the directory.")
        return

    # Loop through each file and process
    for file in metaphlan_files:
        try:
            # Extract sample name from the file path (modify according to your filename format)
            sample_name = os.path.basename(file).replace(".txt", "")
            print(f"Processing file: {file}, Sample name: {sample_name}")
            
            # Read the file
            with open(file, 'r') as f:
                lines = f.readlines()
            
            
            # Load data into a pandas DataFrame
            data = pd.read_csv(file, sep="\t", comment='#',
                                names=['clade_name', 'clade_taxid', 'relative_abundance', 'estimated_number_of_reads_from_the_clade'],
                                header=None,
                                index_col=False
                                )
            
            
            
            # Add sample_name column
            data['sample_name'] = sample_name

            # Append to the list of data frames
            dfs.append(data)
            

        except Exception as e:
            print(f"Error processing file {file}: {e}")

    if not dfs:
        print("No valid data found in any files.")
        return

    # Combine all data frames into a single DataFrame
    combined_df = pd.concat(dfs, ignore_index=True)

    # Reorder columns and save to a TSV file
    combined_df = combined_df[['clade_name', 'clade_taxid', 'relative_abundance', 'estimated_number_of_reads_from_the_clade', 'sample_name']]
    combined_df.to_csv(output_file, sep="\t", index=False)
    print(f"Combined results saved to {output_file}")

if __name__ == "__main__":
    main()