import numpy as np
import sys
import csv

def extract_cwz_and_save(input_file, output_file):
    # Load the data from the .npz file
    data = np.load(input_file, allow_pickle=True, encoding='latin1')
    
    # Extract the 'results_cwz' array
    results_cwz = data['results_cwz']

    # Save the 'results_cwz' array to a CSV file
    with open(output_file, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        # Writing each value in a new row
        for value in results_cwz:
            writer.writerow([value])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python wisecondor_extract_cwz_bowtie2_b37.py <input_file.npz> <output_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_cwz_and_save(input_file, output_file)
