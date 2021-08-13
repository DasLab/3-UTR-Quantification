import pandas as pd

def extract_gene_name(metadata):
    '''function take a section of string in df column called gene_name and splits into a column only what is followed after the label input string in column output is section of string in new column'''
    if len(metadata.split('gene_name')) > 1:
        return metadata.split('gene_name')[1].split('; ')[0].strip().strip('"')
    else:
        return None

def extract_gene_id(metadata):
    if len(metadata.split('gene_id')) > 1:
        return metadata.split('gene_id')[1].split('; ')[0].strip().strip('"')
    else:
        return None
    
def extract_exon_number(metadata):
    if len(metadata.split('exon_number')) > 1:
        return metadata.split('exon_number')[1].split('; ')[0].strip().strip('"')
    else:
        return None
    
    
def extract_gene_biotype(metadata):
    if len(metadata.split('gene_biotype')) > 1:
        return metadata.split('gene_biotype')[1].split('; ')[0].strip().strip('"')
    else:
        return None
    
    
def extract_bed(bedfile_path):
    
    df = pd.read_csv(bedfile_path, sep='\t', header=None) # looking at bed file in system

    df_cols = ['chrom', 'chromstart', 'chromend', 'name', 'score', 'strand', 'chrom', 'source', 'feature', 'thickstart', 'thickend', 'score/strand', 'strand', 'score', 'geneinfo', 'distance']
    df.columns = df_cols
    
    df['gene_name'] = df['geneinfo'].apply(extract_gene_name) # calls on function and looks into column by applying input of what is being parsed
    df['gene_id'] = df['geneinfo'].apply(extract_gene_id)
    df['exon_number'] = df['geneinfo'].apply(extract_exon_number)
    df['gene_biotype'] = df['geneinfo'].apply(extract_gene_biotype)
    
    
    df = df[abs(df.distance) <= 10000] # data looks for column that has values > 10000 and filters into df
    
    return df


# For gene specific dataframe
def count_peaks(df, gene_name):
    fdf = df[df['gene_name'] == gene_name]
    df_hist_y, df_hist_x = np.histogram(abs(fdf['distance']), bins=1000, range=(0, 10000))
    peaks, __ = find_peaks(df_hist_y, distance=150)
    
    filtered_peaks = []
    
    for peak in peaks:
        if df_hist_y[peak] > 0.001*sum(df_hist_y) and df_hist_y[peak] >= 5:
            filtered_peaks.append(peak)
    
    plt.figure(figsize=(20,10))
    plt.plot(df_hist_y)
    plt.plot(peaks, df_hist_y[peaks], 'x')
    plt.plot(filtered_peaks, df_hist_y[filtered_peaks], 'o')
    plt.title(list(df['gene_name'])[0])
    
    return filtered_peaks


# For experiment wide dataframe
def quantify_top_gene(df):
    fdf = df.groupby('gene_name').filter(lambda x: len(x) > 1000)

    unique_gene_names = list(set(fdf['gene_name']))

    for gene_name in unique_gene_names[:10]:
        gene_df = fdf[fdf['gene_name'] == gene_name]

        peaks = count_peaks(gene_df)
        return gene_name, peaks, set(gene_df['exon_number'])
        