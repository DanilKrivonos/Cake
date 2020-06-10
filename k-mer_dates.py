from pyteomics.fasta import read
import pandas as pd
import os
from os.path import split
import argparse
from itertools import permutations

parser = argparse.ArgumentParser(description="Making Spearman correlation")
parser.add_argument('-genome', type=str, help='Genomes file', required=True)
parser.add_argument('-save_way', type=str, help='Output file', required=True)
parser.add_argument('-k', type=int, help='print_kmer', default=2)
args = parser.parse_args()

homology = os.listdir(args.genome)
kmer = args.k
controler = []
nucleotide = ['A', 'G', 'T', 'C']
def mixer(r):
    return set([''.join(i) for i in permutations(nucleotide * r, r=r)])
if kmer > 1:
    controler = []
    for k in mixer(kmer):
        if k in controler:
            continue
        k1 = k
        k2 = k[:: -1]
        if k1[0] == k1[1]:
            if k1[0] == 'A':
                k2 = 'TT'
            elif k1[0] == 'T':
                k2 = 'AA'
            elif k1[0] == 'C':
                k2 = 'GG'
            elif k1[0] == 'G':
                k2 = 'CC'
        controler.append(k1)
        controler.append(k2)
        print('Analizing of {} compounds'.format(k))
        for miorg in homology:
            if '.DS_Store' in miorg:
                continue

            if 'LCB' in miorg:
                print('Analyzing {} ...'.format(miorg[4: -6]))
                save_file = open(args.save_way + '{}_{}&{}.txt'.format(
                    miorg[4: -6], k, k[: :-1]), 'w')
            elif 'SB' in miorg:
                print('Analyzing {} ...'.format(miorg[3: -6]))
                save_file = open(args.save_way + '{}_{}&{}.txt'.format(
                    miorg[3: -6], k, k[: :-1]), 'w')
            else:
                print('Analyzing ...')
                save_file = open(args.save_way, 'w')
            fasta = read(args.genome + miorg)
            for SB in fasta:
                length = len(SB.sequence)
                if length > 1000:
                    rich = ((SB.sequence.count(k[: :-1]) + SB.sequence.count(k)) * 100) / (length - (len(k) - 1))
                    save_file.write(split(split(SB.description)[0])[1] + '\t' + str(len(SB.sequence)) + '\t' + str(rich) + '\n')

            save_file.close()
            print('Text file was saved!')
else:
    for k in mixer(kmer):
        if k in controler:
            continue
        print('Analizing of {} compounds'.format(k))
        for miorg in homology:
            if '.DS_Store' in miorg:
                continue

            if 'LCB' in miorg:
                print('Analyzing {} ...'.format(miorg[4: -6]))
                save_file = open(args.save_way + '{}_{}.txt'.format(
                    miorg[3: -6], k), 'w')
            elif 'SB' in miorg:
                print('Analyzing {} ...'.format(miorg[3: -6]))
                save_file = open(args.save_way + '{}_{}.txt'.format(
                    miorg[2: -6], k), 'w')
            else:
                print('Analyzing ...')
                save_file = open(args.save_way, 'w')
            fasta = read(args.genome + miorg)
            for SB in fasta:
                length = len(SB.sequence)
                if length > 1000:
                    rich = (SB.sequence.count(k) * 100) / (length - (len(k) - 1))
                    save_file.write(split(split(SB.description)[0])[1] + '\t' + str(len(SB.sequence)) + '\t' + str(rich) + '\n')

            save_file.close()
            print('Text file was saved!')

print('Done!')