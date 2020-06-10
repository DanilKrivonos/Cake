import argparse
import pandas as pd
from scipy.stats import wilcoxon
from scipy.stats import kstest
from scipy.stats import mannwhitneyu
from pyteomics.fasta import read
import matplotlib.pyplot as plt
from statistics import median
from os import listdir
import sys
plt.style.use('seaborn')

parser = argparse.ArgumentParser('Mann-whitney test')
parser.add_argument('-all_kmers', 
    help='Path to organisms', 
    required=True
    )
parser.add_argument('-image_save', 
    help='Path to save fold',  
    required=True
    )
args = parser.parse_args()
organisms = listdir(args.all_kmers)
all_2mers = {}
for org in organisms: 
    if '.DS_Store' in org:
        continue
    print(org)
    output_file = pd.read_csv('{}{}/{}_2mers.txt'.format(args.all_kmers, org, org), header=None, names=('Ratio', 'Num'), sep='\t')
    for ind in range(len(output_file)):
        if 'genome' in  output_file['Ratio'][ind]:
            if str(output_file['Ratio'][ind]).split()[-1] not in all_2mers:
                all_2mers[str(output_file['Ratio'][ind]).split()[-1]] = []
                all_2mers[str(output_file['Ratio'][ind]).split()[-1]].append(output_file['Num'][ind + 1])
            else:
                all_2mers[str(output_file['Ratio'][ind]).split()[-1]].append(output_file['Num'][ind + 1])

fig, axes = plt.subplots(nrows=1, figsize=(15,8))
axes.yaxis.grid(True)
axes.set_ylabel('Ratio of frequency 2-mers, end/mid')
boxploter = []
labels = []
sort_median = []
for key in all_2mers:
    labels.append(key)
    boxploter.append(all_2mers[key])
    sort_median.append(median(all_2mers[key]))
sort_median.sort()
sort_boxploter = []
sort_lables = []
for sor_bp in sort_median:
    for ind in range(len(boxploter)):
        if median(boxploter[ind]) == sor_bp:
            sort_boxploter.append(boxploter[ind])
            sort_lables.append(labels[ind])

plt.boxplot(sort_boxploter, labels=sort_lables)

plt.savefig(args.image_save + 'BOX_PLOT_2mers{}.pdf'.format(args.all_kmers.split('/')[-1]), format='pdf', dpi=300)