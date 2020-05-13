import argparse 
import pandas as pd
import numpy as np
import sys
from pyteomics.fasta import read
import os
from os.path import splitext 
from os.path import split

parser = argparse.ArgumentParser(description="Let's make this blocks")
parser.add_argument('-orgs', type=str, help='Prgs', required=True)
parser.add_argument('-coord', type=str, help='Coordinate', required=True)
parser.add_argument('-save_way', type=str, help='Path_of_save', required=True)
parser.add_argument('-syn_fa', type=str, help='Chose your synteny finder', required=True)
parser.add_argument('-out', type=str, help='Output file of synteny finder', default=None, required=True)
parser.add_argument('-refpa', type=str, help='For Sibelia fold with bacterias', default=None)
args = parser.parse_args()

database = os.listdir(args.orgs)
org_count = 0
print('Starting to bake your cake...')
for org in database:
    if '.DS_Store' not in org:

        #list of stains
        if args.syn_fa == 'Mauve':
            if args.out != None:  
                stains = []
                #Change path if you have enother
                path = open(args.out + org + '_results/' + org +'.out')
                for line in path:
                    if line[: 9] != 'terminate':
                        if line[: 3] == '../':
                            stains.append(split(splitext(splitext(line)[0])[0])[1])
                    else:
                        break
        elif args.syn_fa == 'Sibelia':
            #Change path if you have enother
            p = pd.read_csv(args.out + '/{}_results_sibelia/blocks_coords.txt'.format(org), header=None, names=('ind', 'size', 'describtion'), sep='\t')
            check_list = {}
            for i in range(len(p)):
                if str(p['describtion'][i]) != 'nan':
                    if str(p['describtion'][i]) != 'Description':
                        check_list[str(p['describtion'][i])] = {}
                else:
                    break

            orgs = os.listdir(args.refpa + org + '/All/')
            
            for sta in orgs:
                if '.DS_Store' in sta:
                    continue
                op_pa = read(args.refpa + org + '/All/{}'.format(sta))
                for line in op_pa:
                      for num in check_list:
                        if num in line.description:
                            check_list[num] = [sta[: -4]]
            
            stains = []
            for ind in check_list:
                    stains.append(check_list[ind][0])

        #chenge path for your data
        if args.syn_fa == 'Mauve':
            coord = pd.read_csv(
                args.coord + org + '_results/' + org + '.txt', 
                sep=' ',
                header=None,
                names=['sb1', 'sb2', 'genome', 'start1', 'end1', 'start2', 'end2', 'len']
                )
        elif args.syn_fa == 'Sibelia':
             coord = pd.read_csv(
                args.coord + org + '_results_sibelia/' + org + '.txt', 
                sep=' ',
                header=None,
                names=['sb1', 'sb2', 'genome', 'start1', 'end1', 'start2', 'end2', 'len']
                )  

        dataset = os.listdir(args.orgs + org + '/All/')
        for refsta in stains:
            print('+ one strain is ready')
            for stain in dataset:
                if '.DS_Store' not in stain:
                    if refsta == stain[: -4]:
                        coord_g = coord.loc[coord.genome == stain[: -4]]
                        fasta = read(args.orgs + org + '/All/' + stain)
                        for line in fasta:
                            seq = line.sequence

                        edge_cake = ''
                        cake = '' 
                        for i in coord_g.index:
                            seq_i = seq[coord_g.start1[i]:coord_g.end1[i]]
                            #Rewrite condition. It depend on aim
                            #if 0.2 * len(seq_i) < 20: for 20 b.p.
                            #if 0.2 * len(seq_i) < 50: for 50 b.p.
                            #if 0.2 * len(seq_i) < 150: for 150 b.p.
                            #if len(seq_i) < 50: for 10%
                            if 0.2 * len(seq_i) < 150:
                                continue
                            left_edge = seq_i[: 150]
                            right_edge = seq_i[len(seq_i) - 150 :]
                            mid = seq_i[150 : len(seq_i) - 150]
                            cake += ('>sb{}'.format(i))
                            cake += ('\n')
                            cake += (seq_i)
                            cake += ('\n')
                            edge_cake += '>sb{}_right\n'.format(i)
                            edge_cake += right_edge + '\n'
                            edge_cake += '>sb{}_left\n'.format(i) 
                            edge_cake += left_edge + '\n'

                        savefile = open(args.save_way + org + '/' + org + '_' + stain[: -4] + '.fasta', 'w')
                        savefile.write(cake)
                        savefile.close()
                        savefile = open(args.save_way  + org + '/' + org + '_' + stain[: -4] + '_edges.fasta', 'w')
                        savefile.write(edge_cake)
                        savefile.close()

    print('+ one organism is ready (ﾟ▽ﾟ)/')

print('Done! (´｡• ω •｡`)')