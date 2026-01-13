import pandas as pd
import sys
import os

def filter_manifest(manifest_path, outlier_path, output_path, outlier_output_path, merged_output_path):
    # Read the manifest file with explicit column names
    # manifest_df = pd.read_csv(manifest_path, sep='\t', header=None, names=['sample_id', 'type', 'path'])
    manifest_df = pd.read_csv(manifest_path, sep='\t')
    manifest_df.columns = ['sample_id', 'type', 'path']

    # Read the outlier file with explicit column names
    outlier_df = pd.read_csv(outlier_path, sep='\t', header=0, names=['sample_id', 'count', 'outlier'])
    
    # Remove the suffix from the outlier sample_id
    outlier_df['sample_id'] = outlier_df['sample_id'].str.replace('.str_profile', '', regex=False)

    # Ensure the output directory exists
    for path in [output_path, outlier_output_path, merged_output_path]:
        output_dir = os.path.dirname(path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

    # Save the outlier dataframe to a file
    print(f"Saving outliers to: {outlier_output_path}")
    outlier_df.to_csv(outlier_output_path, sep='\t', index=False)

    # Merge the manifest with the outlier information
    merged_df = pd.merge(manifest_df, outlier_df[['sample_id', 'outlier']], on='sample_id', how='left')

    # Save the merged dataframe to a file
    print(f"Saving merged dataframe to: {merged_output_path}")
    merged_df.to_csv(merged_output_path, sep='\t', index=False)

    # Ensure 'outlier' is properly formatted
    merged_df['outlier'] = merged_df['outlier'].astype(str).str.upper()

    # Filter the manifest based on the outlier column
    filtered_df = merged_df[merged_df['outlier'] != 'TRUE']

    # Drop the 'outlier' column
    filtered_df = filtered_df.drop(columns=['outlier'])

    # Save the filtered manifest to the output path without the header
    print(f"Saving filtered manifest to: {output_path}")
    filtered_df.to_csv(output_path, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <outlier_path> <manifest_path> <output_path> <outlier_output_path> <merged_output_path>")
        sys.exit(1)

    outlier_path = sys.argv[1]
    manifest_path = sys.argv[2]
    output_path = sys.argv[3]
    outlier_output_path = sys.argv[4]
    merged_output_path = sys.argv[5]

    filter_manifest(manifest_path, outlier_path, output_path, outlier_output_path, merged_output_path)