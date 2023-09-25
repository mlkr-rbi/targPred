import numpy as np 
import pandas as pd
import openpyxl 
import pyranges as range

'''The original txt file has the headers 
and exon entries removed. This program converts the txt file into a numpy array. 
Entries that do not pass the quality threshold are removed. Then, the pos ID are 
assigned to a range for a locus, which is added as a separate column. 
'''
# import txt file with no header or exons as tsv file 
def load_tab_delimited_file(file_path):
    data = np.genfromtxt(file_path, delimiter='\t', dtype = None, encoding = 'utf-8')
    return data

# Replace 'input.txt' with the path to your tab-delimited TXT file
txt_file_path = '/Users/katarinapejcinovic/Documents/IRB/vcf_files/Klegzom_intron/13HH_S13_no_header_or_exon.txt'

# Load the tab-delimited file into a NumPy array
np_vcf = load_tab_delimited_file(txt_file_path)

df = pd.DataFrame(np_vcf)

# remove SNPs that are poor quality
high_qual = df[df.iloc[:, 6] == "PASS"]

#remove unnecessary columns
clean_vcf = high_qual.drop(high_qual.columns[6:10], axis = 1 )
clean_vcf.columns = ['CHROM',	'POS',	'ID',	'REF',	'ALT', 'QUAL']

# Replace 'output.xlsx' with the desired Excel file path
excel_file_path = 'view_vcf_numpy.xlsx'

# Export the DataFrame to an Excel file
clean_vcf.to_excel(excel_file_path, index = False, header= True)



#This is for making an excel sheet so it's easier to view 
# Convert the NumPy array to a pandas DataFrame
#df = pd.DataFrame(numpy_array)



