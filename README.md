[Purpose] This script processes cell phenotype data to calculate and classify spatial interactions between antigen-presenting cells (APCs) and T cells within regions of interest (ROIs). It identifies contacts such as dyads (APC and CD4T contact) and triads (APC in contact with both CD4T and CD8T), providing a quantitative output for downstream analysis.

[Description]
The script is designed for high-throughput analysis of phenotypic spatial cell data, particularly in immunology contexts. It performs the following steps:
 1. Input Validation: Ensures that the required columns exist in the input data.
 2. Cell Type Generalization: Maps specific T cell phenotypes into broader categories like CD4T, CD8T, or APC.
 3. ROI-based Processing: Splits the dataset by regions of interest (ROIs) and processes each ROI in parallel to improve efficiency.
 4. APC Identification: Identifies antigen-presenting cells based on marker expressions (CD11c and HLADR) and their cell types.
   1) Contact Calculation: Calculates distances between APCs and T cells (CD4T, CD8T) to classify contacts:
   2) Triad: APC in contact with both CD4T and CD8T cells.
   3) Dyad: APC in contact with only CD4T cells.
   4) No Contact: No significant contact.
    Output Generation:
  5. Total results (01.total_results.csv): All analyzed APCs with contact classifications.
  6. Filtered results (02.filtered_results.csv): Only APCs with significant contacts (dyads or triads).


Usage: triad_dyad.py [-h] input_csv [output_csv_total] [output_csv_filtered]

Process phenotype data and calculate cell contacts.

positional arguments:
  input_csv            Input CSV file with raw phenotype data
  output_csv_total     Output CSV file for total results (default: 01.total_results.csv)
  output_csv_filtered  Output CSV file for filtered results (default: 02.filtered_results.csv)

Example:
1. python triad_dyad.py input.csv
2. python triad_dyad.py input.csv output_total output_filter
3. pytnon triad_dyad.py input.csv output_total.csv output_filter.csv
