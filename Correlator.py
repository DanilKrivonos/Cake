import scipy as sp
from scipy.stats import spearmanr
from pyteomics.fasta import read
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import os
from os.path import split
import sys
import numpy as np
import argparse
from itertools import permutations
from sklearn.metrics import r2_score
from scipy.stats import poisson
plt.style.use('seaborn')

parser = argparse.ArgumentParser(description="Making Spearman correlation")
parser.add_argument('-genome', 
    type=str, 
    help='Genomes file', 
    required=True
    )
parser.add_argument('-save_way', 
    type=str, 
    help='Output file', 
    required=True
    )
parser.add_argument('-type', 
    type=str, 
    help='Analyze type', 
    required=True
    )
parser.add_argument('-n1', 
    type=str, 
    help='Nucleotite compaunds', 
    required=True
    )
parser.add_argument('-n2', 
    type=str, 
    help='Nucleotite compaunds', 
    required=True
    )
args = parser.parse_args()

homology = os.listdir(args.genome)
nuc1 = args.n1
nuc2 = args.n2

if args.type == 'full':
    for miorg in homology:
        if '.DS_Store' in miorg:
            continue
        x = []
        y = []
        if 'LCB' in miorg:
            print('Analyzing {} ...'.format(miorg[4: -6]))
            save_file = open(args.save_way + '{}/{}_{}&{}.txt'.format(miorg[4: -6], miorg[4: -6], nuc1, nuc2), 'w')
        elif 'SB' in miorg:
            print('Analyzing {} ...'.format(miorg[3: -6]))
            save_file = open(args.save_way + '{}/{}_{}&{}.txt'.format(miorg[3: -6], miorg[3: -6], nuc1, nuc2), 'w')
        else:
            print('Analyzing ...')
            save_file = open(args.save_way, 'w')
        fasta = read(args.genome + miorg)
        for SB in fasta:
            if len(SB.sequence) > 1000:
                length = len(SB.sequence)
                rich = ((SB.sequence.count(nuc1) + SB.sequence.count(nuc2)) * 100) / (length - (len(nuc1) - 1))
                save_file.write(split(split(SB.description)[0])[1] + '\t' + str(len(SB.sequence)) + '\t' + str(rich) + '\n')
                x.append(length)
                y.append(rich)

        save_file.close()
        print('Text file was saved!')
        plt.figure(figsize=(10, 6))
        plt.scatter(x, y, color='red')
        plt.xlabel("Legth of SB, b.p.", labelpad=13)
        plt.ylabel("{} & {} abundance, %".format(nuc1, nuc2), labelpad=13)
        if 'LCB' in miorg:
            plt.title("{}".format(miorg[4:-6]))
        elif 'SB' in miorg:
            plt.title("{}".format(miorg[3:-6]))
        if 'LCB' in miorg:
            plt.savefig(args.save_way + '{}/{}_{}&{}.pdf'.format(miorg[4:-6], miorg[4:-6], nuc1, nuc2), dpi=300)
        elif 'SB' in miorg:
            plt.savefig(args.save_way + '{}/{}_{}&{}.pdf'.format(miorg[3: -6], miorg[3: -6], nuc1, nuc2), dpi=300)
        else:
            plt.savefig(args.save_way, dpi=300)
else:
    x = []
    y = []
    save_file = open(args.save_way + '{}_{}&{}_every.txt'.format(args.save_way.split('/')[-2], nuc1, nuc2), 'w')
    
    for miorg in homology:
        if '.DS_Store' in miorg:
            continue
        if 'LCB' in miorg:
            print('Analyzing {} ...'.format(miorg[4: -6]))
        elif 'SB' in miorg:
            print('Analyzing {} ...'.format(miorg[3: -6]))
        else:
            print('Analyzing ...')
        if '.DS_Store' in miorg:
            continue
        fasta = read(args.genome + miorg)
        for SB in fasta:
            if 'N' in SB.sequence:
                continue
            if len(SB.sequence) > 1000:
                length = len(SB.sequence)
                rich = ((SB.sequence.count(nuc1) + SB.sequence.count(nuc2)) * 100) / (length - (len(nuc1) - 1))
                save_file.write(split(split(SB.description)[0])[1] + '\t' + str(len(SB.sequence)) + '\t' + str(rich) + '\n')
                x.append(length)
                y.append(rich)

    save_file.close()
    print('Text file was saved!')
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, color='red')
    plt.xlabel("Legth of SB, b.p.", labelpad=13)
    plt.ylabel("{} & {} abundance, %".format(nuc1, nuc2), labelpad=13)
    plt.title("Enterobacteriaceae")
    if 'LCB' in miorg:
        plt.savefig(args.save_way + '{}_{}&{}_every.pdf'.format(args.save_way.split('/')[-2], nuc1, nuc2), dpi=300)
    elif 'SB' in miorg:
        plt.savefig(args.save_way + '{}_{}&{}_every.pdf'.format(args.save_way.split('/')[-2], nuc1, nuc2), dpi=300)
    else:
        plt.savefig(args.save_way, dpi=300)
        
print('Done!')