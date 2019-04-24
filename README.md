
# Upper Quartile Normalization of Raw Count Data From TCGA

### See the `normalization.ipynb` file for the interactive notebook

This guide shows how to reproduce the normalization technique used by The Cancer Genome Atlas (TCGA) to product the data labeled “illuminahiseq_rnaseqv2-RSEM_genes_normalized” from the Brad Institute Firehose (https://gdac.broadinstitute.org/)

I was motivated to make this because I had a hard time googling around to figure out how this is done. This is 
partially due to deprecated links associated with the TCGA. This guide is not perfect and may contain errors. I am open to feedback. 

The “illuminahiseq_rnaseqv2-RSEM_genes_normalized” data from the firehose is produced by performing an upper quartile based normalization on the  preprocessed data in the  “PAAD.uncv2.mRNAseq_raw_counts.txt” file. The way I understand it, upper quartile normalization is used to negate some of the issues associated with a few highly expressed genes “stealing” flow-cell space that would otherwise to reads associated with more lowly expressed genes. There is a good discussion about this here:

https://www.biostars.org/p/314454/

Prerequisites

1. After cloning this repository, you are going to need to navigate to https://gdac.broadinstitute.org/ and click on the browse link under the data heading. I am going to use Pancreatic adenocarcinoma data as an example. 

2. Click on the mRNASeq bar to bring up the mRNASeqArchives

3. Download the mRNAseq_Preprocess  link. Make sure it is the gzipped file and not the MD5 file.

4.  Download the luminahiseq_rnaseqv2-RSEM_genes_normalized link. Make sure it is the gzipped file and not the MD5 file.

4. Extract the file associated with mRNAseq Preprocess link that you downloaded. Understand that “PAAD.uncv2.mRNAseq_raw_counts.txt” should be the exact same file that is already in the data folder of the repository. This is the raw data that we are trying to normalize.

5. Extract the file associated with the  luminahiseq_rnaseqv2-RSEM_genes_normalized link. Understand the “PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt” file contains the normalized data that we are trying to reproduce. Similarly there is a copy of this file already in the data folder of the repository. 


```python
# Import pandas and numpy (you may need to install these first)
import numpy as np
import pandas as pd
```


```python
# Read in the raw counts as a pandas dataframe
not_normalized_df = pd.read_csv('data/PAAD.uncv2.mRNAseq_raw_counts.txt', sep='\t')
```


```python
# Visualize the head of the file to make sure it read in correctly
not_normalized_df.head()
```


```python
# Calculate the 75th percentile for all columns containing raw count data
# The 75th percentile is a number where 75% of the values in each column fall below that number
# More can be read bout percentiles at wikipedia: https://en.wikipedia.org/wiki/Percentile
# We are also going to remove all zero values before calculating the 75th percentile by replacing 0's with NaN, then dropping NaN values
# We are going to save these 75th percentile values in a dictionary
not_normalized_zero_replaced_df = not_normalized_df.replace(0, np.nan)
percentile_dict = {}
for column in not_normalized_df.columns.values[1:]:
    series = not_normalized_zero_replaced_df[column]
    series = series[~series.isnull()]
    percentile = np.percentile(series, 75)
    percentile_dict[column] = percentile
```


```python
# Visualize the percentile dictionary
percentile_dict
```


```python
# Now we normalize the raw count data by dividing each raw count value by the 75th percentile of it's respective column
# we then multiply this quotient by 1000
normalized_data_df = pd.DataFrame()
normalized_data_df['HYBRIDIZATION R'] = not_normalized_df['HYBRIDIZATION R']
for column in not_normalized_df.columns.values[1:]:
    series = not_normalized_df[column]
    normalized_data_df[column] = (series / percentile_dict[column]) * 1000
    
```


```python
# Visualize the data that we normalized
normalized_data_df.head(10)
```


```python
# Read in the normalized counts as a pandas dataframe
normalized_df = pd.read_csv('data/PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep='\t')

# Drop the first row we don't need it
normalized_df = normalized_df[1:]
```


```python
# Visualize the normalized counts as calculated by TCGA
# Comparing this to our normalization, the values are extremely close. There does appear to be a slight difference. 
# My guess is that this has to do with how numpy calculates percentiles, but I am unsure. I would love to hear input from anyone that may know why.
normalized_df.head(10)
```


```python

```
