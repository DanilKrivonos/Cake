import argparse
import numpy as np
import pandas as pd
import csv
from pyteomics.fasta import read
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import sys
plt.style.use('seaborn')

parser = argparse.ArgumentParser('Mann-whitney test')
parser.add_argument('-stain', 
    help='Path to organisms', 
    required=True
    )
parser.add_argument('-coord', 
    help='Path to coordinate', 
    required=True
    )
parser.add_argument('-save_way', 
    help='Path to save fold',  
    required=True
    )
args = parser.parse_args()
print('Starting of analyze ...')

fasta = read(args.stain)
coord = open(args.coord)

read_tsv = pd.DataFrame (csv.reader(coord, delimiter="\t"))
read_tsv = read_tsv.T.drop(0)
Methyla_coordinates = []
for ind in range(len(read_tsv)):
    if ind == 0:
        continue
    Methyla_coordinates.append([int(read_tsv[0][ind]), int(read_tsv[1][ind])])

SB_coordinates = []
SBs = read(args.stain)
SB_dict = {}
for line in SBs:
    if 'GCF_000005845.2_ASM584v2' in line.description:
        SB_dict[line.description.split()[0]] = {}
        SB_dict[line.description.split()[0]]['Coordinates'] = []
        if 'SB' in args.stain:
            SB_dict[line.description.split()[0]]['Coordinates'].append(int(line.description.split()[3].split(':')[0][1:]))
            SB_dict[line.description.split()[0]]['Coordinates'].append(int(line.description.split()[3].split(':')[1][: -1]))
            SB_dict[line.description.split()[0]]['Sequence'] = line.sequence
        else:
            SB_dict[line.description.split()[0]]['Coordinates'].append(int(line.description.split()[2].split(':')[0][1:]))
            SB_dict[line.description.split()[0]]['Coordinates'].append(int(line.description.split()[2].split(':')[1][: -1]))
            SB_dict[line.description.split()[0]]['Sequence'] = line.sequence

print('Analyzing methylome compaunds ...')
Methylation_sites = []
freq_list_mid = []
freq_list_end = []
for SB in SB_dict:
    count_mid = 0
    count_end = 0
    for coord in Methyla_coordinates:
        if SB_dict[SB]['Coordinates'][0] + 150 < coord[0] and SB_dict[SB]['Coordinates'][1] - 150 > coord[1]:
            count_mid += 1
        elif SB_dict[SB]['Coordinates'][0] < coord[0] and SB_dict[SB]['Coordinates'][0] + 150 > coord[1]:
            count_end += 1
        elif SB_dict[SB]['Coordinates'][1] - 150 < coord[0] and SB_dict[SB]['Coordinates'][1] > coord[1]:
            count_end += 1
    if count_mid != 0:
        freq_list_mid.append(count_mid)
    if count_end != 0:
        freq_list_end.append(count_end)

fig, axes = plt.subplots(nrows=1, figsize=(15,8))
axes.yaxis.grid(True)
axes.set_ylabel('Frequency of methylated sites, %')
savefile = open(args.save_way + 'RESULT_METHYLOME_SIBELIA.txt', 'w')
savefile.write('Mann Whitneyu results\tP-value = {}\n'.format(mannwhitneyu(freq_list_end, freq_list_mid)[1]))
savefile.close()
print('P-value = {}'.format(mannwhitneyu(freq_list_end, freq_list_mid)[1]))
plt.boxplot([freq_list_end, freq_list_mid], labels=['END', 'MID'])
plt.savefig(args.save_way + 'BOX_PLOT_METHYLOME_SIBELIA.pdf', format='pdf', dpi=300)
print('Done!')