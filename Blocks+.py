import argparse 
import pandas as pd
import numpy as np
import sys
from pyteomics.fasta import read

parser = argparse.ArgumentParser(description="Let's make this blocks")
parser.add_argument('-fasta', type=str, help='fasta', required=True)
parser.add_argument('-coord', type=str, help='coordinate', required=True)
parser.add_argument('-save_way', type=str, help='path_of_save', required=True)
args = parser.parse_args()

fasta = read(args.fasta)
coord = pd.read_csv(args.coord, sep=' ', header=None, names=['sb1', 'sb2', 'genome', 'start1', 'end1', 'start2', 'end2', 'len'])

coord_g1 = coord.loc[coord.genome == 'g1']

for line in fasta:
    seq = line.sequence
edge_cake = ''
cake = '' 
for i in coord_g1.index:
    seq_i = seq[coord_g1.start1[i]:coord_g1.end1[i]]
    if len(seq_i) < 102:
        continue
    left_edge = seq_i[: 150]
    right_edge = seq_i[len(seq_i) - 150 :]
    mid = seq_i[150 : len(seq_i) - 150 ]
    cake += ('>sb{}'.format(i))
    cake += ('\n')
    cake += (seq_i)
    cake += ('\n')
    edge_cake += '>sb{}_right\n'.format(i)
    edge_cake += right_edge + '\n'
    edge_cake += '>sb{}_left\n'.format(i) 
    edge_cake += left_edge + '\n'

savefile = open(args.save_way + '.fasta', 'w')
savefile.write(cake)
savefile.close()

savefile = open(args.save_way + 'edges.fasta', 'w')
savefile.write(edge_cake)
savefile.close()