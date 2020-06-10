from itertools import permutations
import argparse 
import sys
import os
from pyteomics.fasta import read
from scipy.stats import chisquare
import math 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy import stats
import matplotlib.pyplot as mpl
from statsmodels.stats.multitest import multipletests
plt.style.use('seaborn')
plt.rcParams["font.family"] = "Times New Roman"

parser = argparse.ArgumentParser(description='Use to find kmers')
parser.add_argument('-modifa', 
    type=str, 
    help='modifying_fasta', 
    required=True
    )
parser.add_argument('-reffa', 
    type=str, 
    help='reference_fasta', 
    required=True
    )
parser.add_argument('-save_way', 
    type=str, 
    help='path_of_save', 
    required=True
    )
parser.add_argument('-image_save', 
    type=str, 
    help='image_save_path', 
    required=True
    )
args = parser.parse_args()
ref_path = args.reffa
        
mpl.rcParams['figure.figsize'] = [12, 8]
print('Analizing of motive abundance of E.coli ...\n')
genome = read(args.reffa)
syntheny_bloks = read(args.modifa)
alphablet = ['A', 'T', 'G' ,'C']

new_genome = '' 
for line in genome:
    new_genome += line.sequence

def create_lib(r):
    return set([''.join(i) for i in permutations(alphablet * r, r=r)])

all_ends = ''
for line in syntheny_bloks:
    all_ends += line.sequence

print('Waiting...\n')

lable_k = []
all_sb_friq = []
end_sb_friq = []
sort_sb_friq = []
sort_end_sb_friq = []
all_p_values = []
sort_p_values = []
sort_lable_k = []
cake = ''

dict_of_motiv = {'GATC': ['GATC'], 'CANCATC': [],'AACN4CTTT': [], 'RTACN4GTG': [], 'GAGACC': ['GAGACC']}
for kmer in create_lib(1):
    dict_of_motiv['CANCATC'].append('CA{}CATC'.format(kmer))
for kmer in create_lib(4):
    dict_of_motiv['AACN4CTTT'].append('AAC{}CTTT'.format(kmer))
    dict_of_motiv['RTACN4GTG'].append('ATAC{}GTG'.format(kmer))
    dict_of_motiv['RTACN4GTG'].append('GTAC{}GTG'.format(kmer))
    

for key in dict_of_motiv:
    lable_k.append(key)
    mid_count = 0
    end_count = 0
    print('Counting {} ...'.format(key))
    for motiv in dict_of_motiv[key]:
        mid_count += new_genome.count(motiv)
        end_count += all_ends.count(motiv)
        len_mot = len(motiv)
    mid = mid_count * 100/(len(new_genome) - len_mot + 1)
    end = end_count * 100/(len(all_ends) - len_mot + 1)
    all_sb_friq.append(mid_count * 100/(len(new_genome) - len_mot + 1))
    end_sb_friq.append(end_count * 100/(len(all_ends) - len_mot + 1))
    cake += 'Frequency in ends {} \t{}\n'.format(motiv, end_sb_friq)
    cake += 'Frequency in genome {} \t{}\n'.format(motiv, all_sb_friq)
    cake += 'Ratio \t{}\n'.format(end/mid)
    pvalue = stats.chi2_contingency([
        [len(all_ends), len(new_genome)],
        [end_count, mid_count]
        ])[1]

    print('P-value = {}\n'.format(pvalue))
    cake += 'P-value for {} \t{}\n'.format(motiv, str(math.log1p(pvalue)))
    all_p_values.append(pvalue)
   
#correction for multiple comparisons
all_p_values = multipletests(all_p_values, method='fdr_bh')[1]
all_p_values = list(np.array(all_p_values))
sort_p_values = list(np.sort(np.array(all_p_values)))
for sornum in sort_p_values:
    for num in range(len(all_p_values)):
        if sornum == all_p_values[num]:
            if all_sb_friq[num] not in sort_sb_friq:
                sort_sb_friq.append(all_sb_friq[num])
                sort_end_sb_friq.append(end_sb_friq[num])
                sort_lable_k.append(lable_k[num])
            
                        
savefile = open(args.save_way, 'w')
savefile.write(cake)
savefile.close()

print('Text file is complete!\n༼ つ ◕_◕ ༽つ\n')
print('Now please wait a picture\n(ﾉ◕ヮ◕)ﾉ*:･ﾟ✧\n')

x = np.arange(len(sort_lable_k))
width = 0.35
fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, sort_sb_friq, width, label='K-mer in the middle of SB, %')
rects2 = ax.bar(x + width/2, sort_end_sb_friq, width, label='K-mer in the ends of SB, %')

ax.set_ylabel('Frequency, %')
ax.set_xlabel('K-mers')
ax.set_title('Escherichia coli')
ax.set_xticks(x)
ax.set_xticklabels(sort_lable_k, rotation='vertical')
ax.legend(loc='center', bbox_to_anchor=(1.15, 0.5))

def autolabel(rects):
    for rect in range(len(rects)):
        if rects1[rect].get_height() > rects2[rect].get_height():
            height = rects1[rect].get_height()
            ax.annotate('log(P-value)=\n{}'.format(round(np.log(sort_p_values[rect]), 3)),
                xy=(rects1[rect].get_x() + rects1[rect].get_width() / 2, height),
                xytext=(0, 1), 
                textcoords="offset points",
                ha='center', va='bottom')
        else:
            height = rects2[rect].get_height()
            ax.annotate('log(P-value)=\n{}'.format(round(np.log(sort_p_values[rect]), 3)),
                xy=(rects1[rect].get_x() + rects1[rect].get_width() / 2, height),
                xytext=(0, 1), 
                textcoords="offset points",
                ha='center', va='bottom')

autolabel(rects1)
fig.tight_layout()
plt.savefig(args.image_save, format='pdf', dpi=300)

print('Analyze of Escherichia coli was done!\n')