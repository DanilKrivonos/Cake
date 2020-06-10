import argparse
from os import listdir
from pyteomics.fasta import read
from itertools import permutations
import sys

parser = argparse.ArgumentParser('Mann-whitney test')
parser.add_argument('-orgs', 
    help='Path to organisms', 
    required=True
    )
parser.add_argument('-save_way', 
    help='Path to save fold',  
    required=True
    )
args = parser.parse_args()

print('Starting of analyze ...')
organisms = listdir(args.orgs)
alphablet = ['A', 'T', 'G' ,'C']
def create_lib(r):
    return set([''.join(i) for i in permutations(alphablet * r, r=r)])

savefile = open(args.save_way + 'FULL_kmers_MAUVE.txt', 'w')
savefile_GATC = open(args.save_way + 'GATC_MAUVE.txt', 'w')
lable_k = []
mid_sb_friq = []
end_sb_friq = []
end_sb_friq_GATC = []
mid_sb_friq_GATC = []
ratio_all_kmers = []
ratio_all_GATC = []
all_p_values_kmers = []
all_p_values_GATC = []

for kmer in create_lib(4):
    if kmer != 'GATC':
        for org in organisms:
            cake = ''
            end_genome = ''
            mid_genome = ''
            if '.DS_Store' in org:
                continue

            print('Analyzin of {} ...'.format(org))
            mid_genome = read('{}{}/{}_middle_SB.fasta'.format(args.orgs, org, org))
            for line in mid_genome:
                mid_genome = line.sequence

            end_genome = read('{}{}/{}_ENDS_SB.fasta'.format(args.orgs, org, org))
            for line in end_genome:
                end_genome = line.sequence

            
            print('Counting {} ...'.format(kmer))
            lable_k.append(kmer)
            mid = mid_genome.count(kmer) * 100/(len(mid_genome) - 5)
            end = end_genome.count(kmer) * 100/(len(end_genome) - 5)
            ratio = end/mid
            ratio_all_kmers.append(ratio)
            cake += 'Frequency in the ends for every 4-mers\t{}\n'.format(end)
            cake += 'Frequency in the middle genome for every 4-mers\t{}\n'.format(mid)
            cake += 'Ratio \t{}\n'.format(ratio)
            savefile.write(cake)

    elif kmer == 'GATC': 
        for org in organisms:
            cake_GATC = ''
            end_genome = ''
            mid_genome = ''
            if '.DS_Store' in org:
                continue

            print('Analyzin of {} ...'.format(org))
            mid_genome = read('{}{}/{}_middle_SB.fasta'.format(args.orgs, org, org))
            for line in mid_genome:
                mid_genome = line.sequence

            end_genome = read('{}{}/{}_ENDS_SB.fasta'.format(args.orgs, org, org))
            for line in end_genome:
                end_genome = line.sequence

            print('Counting {} ...'.format(kmer))
            lable_k.append(kmer)
            mid_GATC = mid_genome.count(kmer) * 100/(len(mid_genome) - 5)
            end_GATC = end_genome.count(kmer) * 100/(len(end_genome) - 5)
            ratio = end_GATC/mid_GATC
            ratio_all_GATC.append(ratio)
            cake_GATC += 'Frequency in the ends for GATC\t{}\n'.format(end_GATC)
            cake_GATC += 'Frequency in the middle genome GATC\t{}\n'.format(mid_GATC)
            cake_GATC += 'Ratio_{}\t{}\n'.format(org, end_GATC/mid_GATC)
            savefile_GATC.write(cake_GATC)

savefile.close()
savefile_GATC.close()