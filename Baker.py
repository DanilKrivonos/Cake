import argparse
import os 
from pyteomics.fasta import read
import sys

parser = argparse.ArgumentParser(description='Use to stick blocks')
parser.add_argument('-fa', type=str, help='Modifying fasta', required=True)
parser.add_argument('-save_way', type=str, help='Path of save', required=True)
args = parser.parse_args()

path_to_org = os.listdir(args.fa)
for org in path_to_org:
    if '.DS_Store' in org:
        continue
    full_dataset = os.listdir(args.fa + org)
    All_SB = ''
    for strain in full_dataset:
        if '.DS_Store' in strain:
            continue
        if 'edges.fasta' not in strain:
            full_genome = read(args.fa + org + '/{}'.format(strain))
            for string in full_genome:
                All_SB += string.sequence
    print('All synteny blocks were combined for {}!'.format(org))
    savefile = open(args.save_way + org + '/{}_middle_SB.fasta'.format(org), 'w')
    savefile.write('>dataset_of_middle_SB')
    savefile.write('\n')
    savefile.write(All_SB)
    savefile.close()
    print('All middle of synteny blocks for {}'.format(org))

for org in path_to_org:
    if '.DS_Store' in org:
        continue
    ends_dataset = os.listdir(args.fa + org)
    all_ends = ''
    for strain in ends_dataset:
        if '.DS_Store' in strain:
            continue
        if 'edges.fasta' in strain:
            ends_gen = read(args.fa + org + '/{}'.format(strain))
            for string in ends_gen:
                all_ends += string.sequence
    print('All ends of synteny blocks were combined {}!'.format(org))
    savefile = open(args.save_way + org + '/{}_ENDS_SB.fasta'.format(org), 'w')
    savefile.write('>dataset_of_ends_SB_ends')
    savefile.write('\n')
    savefile.write(all_ends)
    savefile.close()
    print('All ends saved for {}'.format(org))

print('Done!:)')