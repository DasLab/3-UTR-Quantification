import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


from extract_bed import extract_bed # fix for py file

def find_gene_specific_histo_peaks(gdf, gene_name, bins=1000, max_length=10000):
    
    df_hist_y, df_hist_x = np.histogram(abs(gdf['distance']), bins=bins, range=(0, max_length))
    
    peaks, __ = find_peaks(df_hist_y, distance=150)
    
    filtered_peaks = []
    
    for peak in peaks:
        if df_hist_y[peak] > 0.001*sum(df_hist_y) and df_hist_y[peak] >= 5:
            filtered_peaks.append(peak)
            
    df_hist_y = df_hist_y/sum(df_hist_y)
    
    rows = []
    for peak in filtered_peaks:
        # TODO: sum under curve for isoform abundance
        #https://stackoverflow.com/questions/7958956/remove-data-points-below-a-curve-with-python 
        rows.append({'gene_name': gene_name, 'isoform_length': peak*(max_length/bins), 'isoform_abundance': df_hist_y[peak]})

    
    return pd.DataFrame(rows)


dfs = []
for gene in unique_gene_names:
    tmp_df = fdf[fdf['gene_name'] == gene]
    dfs.append(find_gene_specific_histo_peaks(tmp_df, gene))

peaks_df = pd.concat(dfs)

peaks_df[peaks_df['gene_name'] == 'Plp1']

