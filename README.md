Process phenotype data and calculate cell contacts.

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
