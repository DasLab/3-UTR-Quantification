"""
This is the main APA analyis pipeline file
"""
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from extract_bed import extract_bed
from extract_peaks import find_gene_specific_histo_peaks

def main():

	if len(sys.argv) < 2:
		print('Please enter BED file output from bedtools closest command!')
		print('USAGE: python apa_pipeline.py /path/to/bedfile.bed')
		quit()

	bed_path = sys.argv[1]
	closest_df = extract_bed(bed_path)

	fdf = closest_df.groupby('gene_name').filter(lambda x: len(x) > 500)
	unique_gene_names = list(set(fdf['gene_name']))


if __name__ == "__main__":
	main()