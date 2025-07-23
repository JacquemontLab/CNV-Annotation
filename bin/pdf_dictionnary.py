#!/usr/bin/env python

# Florian Bénitière 16/03/2025
# Script to generate a .pdf that describe the content by column of a .parquet


import sys

import duckdb
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import os
from math import floor  # For rounding down numbers

def generate_pdf_dictionary_duckdb(path_to_snv_dataset, total_memory):
    output_file = path_to_snv_dataset.replace(".parquet", "_dictionary.pdf")

    if os.path.isfile(path_to_snv_dataset) and path_to_snv_dataset.endswith(".parquet"):
        parquet_pattern = path_to_snv_dataset
    elif os.path.isdir(path_to_snv_dataset):
        # Look for subdirectories
        has_subdirs = any(os.path.isdir(os.path.join(path_to_snv_dataset, entry))
                          for entry in os.listdir(path_to_snv_dataset))
        if has_subdirs:
            parquet_pattern = os.path.join(path_to_snv_dataset, "**/*.parquet")
        else:
            parquet_pattern = os.path.join(path_to_snv_dataset, "*.parquet")
    else:
        raise ValueError(f"Input path '{path_to_snv_dataset}' is neither a .parquet file nor a directory.")


    # Open DuckDB connection
    con = duckdb.connect()
    con.execute(f"SET memory_limit = '{total_memory}GB'")
    
    # Get column names and types
    #schema = con.execute(f"DESCRIBE FROM read_parquet('{path_to_snv_dataset}/*/*.parquet')").fetchall()
    schema = con.execute(f"DESCRIBE FROM read_parquet('{parquet_pattern}')").fetchall()


    with PdfPages(output_file) as pdf:
        with tqdm(total=len(schema)) as pbar:
            for column_name, duck_type, *_ in schema:
                print(column_name)
                try:
                    if duck_type in ['INTEGER', 'DOUBLE', 'BIGINT', 'REAL', 'FLOAT']:
                        # Min/max for label
                        min_val, max_val = con.execute(f"SELECT MIN({column_name}), MAX({column_name}) FROM read_parquet('{parquet_pattern}')").fetchone()
                        
                        bin_count = 50
                        bin_width = (max_val - min_val) / bin_count
                        if bin_width == 0:
                            bin_width = 1
                            
                        # Manually compute bin indices
                        histogram_query = f"""
                        SELECT 
                            FLOOR(({column_name} - {min_val}) / {bin_width}) AS bin_idx,
                            COUNT(*) AS count
                        FROM read_parquet('{parquet_pattern}')
                        WHERE {column_name} IS NOT NULL
                        GROUP BY bin_idx
                        ORDER BY bin_idx
                        """
                        bins_data = con.execute(histogram_query).fetchall()
                        if not bins_data:
                            pbar.update(1)
                            continue

                        # Convert bucket numbers to midpoints
                        bin_edges = [min_val + (i * (max_val - min_val) / 50) for i in range(51)]
                        bucket_ids, counts = zip(*bins_data)
                        bin_centers = [bin_edges[int(b)] for b in bucket_ids]

                        # Plot histogram
                        plt.bar(bin_centers, counts, width=(bin_width), color='lightblue', edgecolor='black')
                        plt.title(f"Distribution Column: {column_name}\nDataType: {duck_type}, Min: {min_val:,}, Max: {max_val:,}",
                                  fontsize=12)
                        plt.xlabel(column_name)
                        plt.ylabel("Frequency")
                        plt.tight_layout()
                        pdf.savefig()
                        plt.close()

                    elif duck_type in ['VARCHAR', 'BOOLEAN','VARCHAR[]']:
                        # Frequency table for top 20 values
                        freq_query = f"""
                        SELECT {column_name}, COUNT(*) AS count
                        FROM read_parquet('{parquet_pattern}')
                        GROUP BY {column_name}
                        ORDER BY count DESC
                        LIMIT 20
                        """
                        freq_data = con.execute(freq_query).fetchall()
                        if not freq_data:
                            pbar.update(1)
                            continue
                        
                        df_counts = pd.DataFrame(freq_data, columns=[column_name, 'count'])
                        
                        # Create a new figure for the table
                        fig, ax = plt.subplots(figsize=(10, 6))  # Increase the figure size for the table
                        fig.suptitle(f"First 20 Occurrences of Column: {column_name}, DataType: {duck_type}",
                                    fontsize=12)
                        
                        # Hide axes
                        ax.axis('off')

                        # Create the table and add it to the figure
                        table = ax.table(cellText=df_counts.values, colLabels=df_counts.columns, loc='center', cellLoc='center')

                        # Customize table style (optional)
                        table.auto_set_font_size(False)
                        table.set_fontsize(8)
                        table.scale(1, 1)  # Scale the table to fit more rows
                        
                        # Use tight_layout to ensure everything fits
                        plt.tight_layout()
                        
                        # Bold the first row (column headers)
                        for (i, j), cell in table.get_celld().items():
                            if i == 0:  # Row 0 is the header row
                                cell.set_text_props(fontweight='bold',fontsize=11)
                        
                        # Save the table figure to the PDF
                        pdf.savefig(fig)
                        plt.close(fig)  # Close the table figure to free memory
                    # Update progress bar
                    else:
                        print(f"Not generating for column {column_name}: {e}")
                except Exception as e:
                    print(f"Error while generating for column {column_name}: {e}")
                pbar.update(1)


if __name__ == "__main__":
    parquet_input = sys.argv[1]
    cpus = int(sys.argv[2])
    mem_per_cpu = float(sys.argv[3])
    total_memory = int(cpus * mem_per_cpu)
    
    print(f"[INFO] Input Parquet: {parquet_input}")
    print(f"[INFO] CPUs: {cpus}, Memory per CPU: {mem_per_cpu} GB")
    print(f"[INFO] Total memory allocated: {total_memory} GB")

    generate_pdf_dictionary_duckdb(parquet_input, total_memory)
