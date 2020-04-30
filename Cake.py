from itertools import permutations
import argparse 
import sys
from pyteomics.fasta import read
import scipy.stats as stats
import math 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn')

parser = argparse.ArgumentParser(description='Use to find kmers')
parser.add_argument('-modifa', type=str, help='modifying_fasta', required=True)
parser.add_argument('-reffa', type=str, help='reference_fasta',)
parser.add_argument('-k', type=int, help='print_kmer', default=2)
parser.add_argument('-save_way', type=str, help='path_of_save', required=True)
parser.add_argument('-image_save', type=str, help='image_save_path', required=True)
args = parser.parse_args()

genome = read(args.reffa)
syntheny_bloks = read(args.modifa)
k = args.k
alphablet = ['A', 'T', 'G' ,'C']

new_genome = '' 
for line in genome:
    new_genome += line.sequence

def create_lib(r):
    return set([''.join(i) for i in permutations(alphablet * r, r=r)])

all_ends = ''
for line in syntheny_bloks:
    all_ends += line.sequence

print('Waiting...')
lable_k = []
all_sb_friq = []
end_sb_fiq = []
all_p_values = []
cake = ''
for kmer in create_lib(k):
    print('Counting {} ...'.format(kmer))
    lable_k.append(kmer)
    all_sb_friq.append(new_genome.count(kmer) * 100/(len(new_genome) - k + 1))
    end_sb_fiq.append(all_ends.count(kmer) * 100/(len(all_ends) - k + 1))
    cake += 'Frequency in ends {} \t{}\n'.format(kmer, str(all_ends.count(kmer) * 100/(len(all_ends) - k + 1)))
    cake += 'Frequency in genome {} \t{}\n'.format(kmer, str(new_genome.count(kmer) * 100/(len(new_genome) - k + 1)))
    cake += 'Ratio \t{}\n'.format(str((all_ends.count(kmer) * 100/(len(all_ends) - k + 1))/((new_genome.count(kmer) * 100)/(len(new_genome) - k + 1))))
    Stat_table = [[len(all_ends), len(new_genome)], [all_ends.count(kmer), new_genome.count(kmer)]]
    oddsratio, pvalue = stats.fisher_exact(Stat_table)
    cake += 'P-value for {} \t{}\n'.format(kmer, str(math.log1p(pvalue)))
    all_p_values.append(pvalue)

savefile = open(args.save_way, 'w')
savefile.write(cake)
savefile.close()
print('Text file is complete!\n༼ つ ◕_◕ ༽つ')
print('Now please wait a picture\n(ﾉ◕ヮ◕)ﾉ*:･ﾟ✧')
x = np.arange(len(lable_k))
width = 0.35
fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, all_sb_friq, width, label='K-mer in all SB, %')
rects2 = ax.bar(x + width/2, end_sb_fiq, width, label='K-mer Ends, %')

ax.set_ylabel('Scores')
ax.set_title('Distribution by K-mers composition')
ax.set_xticks(x)
ax.set_xticklabels(lable_k)
ax.legend()

def autolabel(rects):
    for rect in range(len(rects1)):
        height = rects[rect].get_height()
        ax.annotate('log P-value =\n{}'.format(round(np.log(all_p_values[rect]), 3)),
            xy=(rects1[rect].get_x() + rects1[rect].get_width() / 2, height),
            xytext=(0, 5), 
            textcoords="offset points",
            ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)
fig.tight_layout()
plt.savefig(args.image_save, format='pdf', dpi=300)
plt.show()