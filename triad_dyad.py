import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
import multiprocessing
import argparse
import time
from tqdm import tqdm

def check_required_columns(data, required_columns):
    """Check for missing required columns in the input DataFrame."""
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

def classify_cell_type(cell_type):
    """Generalize cell types into broader categories."""
    if 'CD4T_' in cell_type:
        return 'CD4T'
    elif 'CD8T_' in cell_type:
        return 'CD8T'
    elif cell_type not in ['CD4T', 'CD8T']:
        return 'APC' if cell_type == 'APC' else cell_type
    return cell_type

def calculate_contacts(apc_row, target_cells):
    """Calculate distances between an APC row and target cells."""
    distances = np.sqrt((target_cells['X_centroid'] - apc_row['X_centroid'])**2 + 
                        (target_cells['Y_centroid'] - apc_row['Y_centroid'])**2)
    return distances <= apc_row['major_axis_length']

def process_roi(test_roi):
    """Process each ROI to identify APC contacts."""
    results = []
    print(f"Processing ROI: {test_roi['new_roi'].iloc[0]}")

    treg = test_roi[test_roi['cell_type'] == 'Treg']
    cd4t = test_roi[test_roi['cell_type'] == 'CD4T']
    cd8t = test_roi[test_roi['cell_type'] == 'CD8T']

    immune_cell_types = ['DC', 'MDSC', 'Other Myeloid', 'B cells', 'NK_cells']

    def is_apc(row):
        """APC identification logic."""
        if (row['CD11c'] > 0 and row['HLADR'] > 0) and (row['cell_type'] in immune_cell_types):
            return True
        return False

    apc = test_roi[test_roi.apply(is_apc, axis=1)]

    for _, apc_row in apc.iterrows():
        contact_cd4t = any(calculate_contacts(apc_row, cd4t))
        contact_cd8t = any(calculate_contacts(apc_row, cd8t))

        if contact_cd4t and contact_cd8t:
            contact_type = 'triad'
        elif contact_cd4t:
            contact_type = 'dyad'
        else:
            contact_type = 'no contact'

        results.append({'APC_ID': apc_row.name, 
                        'Cell_Type': apc_row['cell_type'], 
                        'ROI': apc_row['new_roi'],
                        'Contact_Type': contact_type})

    print(f"ROI {test_roi['new_roi'].iloc[0]} processed.")
    return pd.DataFrame(results, columns=['APC_ID', 'Cell_Type', 'ROI', 'Contact_Type'])

def ensure_csv_extension(filename):
    """Ensure the file name has a .csv extension."""
    return filename if filename.endswith('.csv') else f"{filename}.csv"

def main(input_csv, output_csv_total=None, output_csv_filtered=None):
    start_time = time.time()

    # Set default output file names if not provided
    if output_csv_total is None:
        output_csv_total = "01.total_results.csv"
    if output_csv_filtered is None:
        output_csv_filtered = "02.filtered_results.csv"

    # Ensure the output file names have .csv extension
    output_csv_total = ensure_csv_extension(output_csv_total)
    output_csv_filtered = ensure_csv_extension(output_csv_filtered)

    # Check if input file exists
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"File '{input_csv}' not found.")

    # Load raw data
    print("Loading data...")
    try:
        data = pd.read_csv(input_csv, low_memory=False, index_col=0)
        print(f"Data loaded: {len(data)} rows")
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    print("First few rows of the loaded data:")
    print(data.head())

    # Generalize cell types
    print("Generalizing cell types...")
    data['cell_type'] = data['cell_type'].apply(classify_cell_type)
    print("Cell type generalization complete.")

    # Preprocess data
    print("Preprocessing data...")
    for col, dtype in {'CD11c': 'int16', 'HLADR': 'int16', 
                       'X_centroid': 'float32', 'Y_centroid': 'float32', 
                       'major_axis_length': 'float32'}.items():
        data[col] = data[col].fillna(0).astype(dtype)
    print("Preprocessing complete.")

    # Process each ROI in parallel
    print("Splitting and processing ROIs...")
    roi_list = data['new_roi'].unique()
    split_data = [data[data['new_roi'] == roi] for roi in roi_list]
    num_cores = min(len(roi_list), multiprocessing.cpu_count())

    print(f"Using {num_cores} cores for parallel processing...")
    results = Parallel(n_jobs=num_cores)(
        delayed(process_roi)(roi_data) for roi_data in tqdm(split_data, desc="Processing ROIs"))

    # Combine results
    print("Combining results...")
    final_results_df = pd.concat(results, ignore_index=True)

    print(f"Total processed rows: {len(final_results_df)}")
    print(final_results_df.head())

    # Save results
    print("Saving results...")
    summary1 = final_results_df
    summary2 = final_results_df[final_results_df['Contact_Type'] != 'no contact']

    summary1.to_csv(output_csv_total, index=False)
    summary2.to_csv(output_csv_filtered, index=False)
    print(f"Results saved to {output_csv_total} and {output_csv_filtered}")

    print(f"Processing completed in {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process phenotype data and calculate cell contacts.')
    parser.add_argument('input_csv', type=str, help='Input CSV file with raw phenotype data')
    parser.add_argument('output_csv_total', type=str, nargs='?', default="01.total_results.csv",
                        help='Output CSV file for total results (default: 01.total_results.csv)')
    parser.add_argument('output_csv_filtered', type=str, nargs='?', default="02.filtered_results.csv",
                        help='Output CSV file for filtered results (default: 02.filtered_results.csv)')
    args = parser.parse_args()

    main(args.input_csv, args.output_csv_total, args.output_csv_filtered)

