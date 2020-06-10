import argparse
import pandas as pd
from scipy.stats import wilcoxon
from scipy.stats import kstest
from scipy.stats import mannwhitneyu
from pyteomics.fasta import read
import matplotlib.pyplot as plt
import sys
plt.style.use('seaborn')

parser = argparse.ArgumentParser('Mann-whitney test')
parser.add_argument('-all_kmers', 
    help='Path to organisms', 
    required=True
    )
parser.add_argument('-GATC', 
    help='Path to organisms', 
    required=True
    )
parser.add_argument('-save_way', 
    help='Path to save fold',  
    required=True
    )
parser.add_argument('-image_save', 
    help='Path to save fold',  
    required=True
    )
args = parser.parse_args()

all_kmers = pd.read_csv(args.all_kmers, header=None,
    names=('Ratio', 'Num'), 
    sep='\t'
    )
GATC =  pd.read_csv(args.GATC, header=None, 
    names=('Ratio', 'Num'), 
    sep='\t'
    )

KMERS_LIST = []
for ind in range(len(all_kmers)):
    if 'Ratio' in all_kmers['Ratio'][ind]:
        KMERS_LIST.append(all_kmers['Num'][ind])

GATC_LIST = []
for ind in range(len(GATC)):
    if 'Ratio' in GATC['Ratio'][ind]:
        GATC_LIST.append(GATC['Num'][ind])

fig, axes = plt.subplots(nrows=1, figsize=(15,8))
axes.yaxis.grid(True)
axes.set_ylabel('Ratio of frequency 4-mers, end/mid')
savefile = open(args.save_way, 'w')
savefile.write('Mann Whitneyu results\tP-value = {}\n'.format(mannwhitneyu(KMERS_LIST, GATC_LIST)[1]))
savefile.close()
print('P-value = {}'.format(mannwhitneyu(KMERS_LIST, GATC_LIST)[1]))
plt.boxplot([KMERS_LIST, GATC_LIST], labels=['All tetramers without GATC', 'GATC'])
plt.savefig(args.image_save, format='pdf', dpi=300)