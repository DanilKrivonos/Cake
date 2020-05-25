import argparse 
import pandas as pd
import numpy as np
import sys
from pyteomics.fasta import read
import os
from os.path import splitext 
from os.path import split

parser = argparse.ArgumentParser(description="Let's make this blocks")
parser.add_argument('-orgs', type=str, help='Orgsnisms paths', required=True)
parser.add_argument('-coord', type=str, help='Coordinate', required=True)
parser.add_argument('-cut', type=str, help='Cut of coordinate', required=True, default=20)
parser.add_argument('-sb_out', type=str, help='If you need SBs fastas',  default='NO')
parser.add_argument('-save_way', type=str, help='Path_of_save', required=True)
parser.add_argument('-syn_fa', type=str, help='Chose your synteny finder', required=True)
parser.add_argument('-out', type=str, help='Output file of synteny finder', required=True)
args = parser.parse_args()

database = os.listdir(args.orgs)
org_count = 0
print('Starting to bake your cake...')
                            
if args.syn_fa == 'Mauve':
    for org in database:
        if '.DS_Store' in org:
            continue
        saveSB = open(args.sb_out + 'LCB_{}.fasta'.format(org), 'w')
        print('Analyzing of {} ...'.format(org))
        #list of stains
        stains = []
        #Change path if you have enother
        path = open(args.out + org + '_results/' + org +'.out')
        for line in path:
            if line[: 9] != 'terminate':
                if line[: 3] == '../':
                    stains.append(split(splitext(splitext(line)[0])[0])[1][: -8])
            else:
                break
        coord = pd.read_csv(
            args.coord + org + '_results/' + org + '.txt', 
            sep='\t',
            header=None,
            names=['sb1', 'sb2', 'genome', 'start', 'end', 'maxlen']
            )

        dataset = os.listdir(args.orgs + org + '/All/')
        for refsta in stains:
            for stain in dataset:
                if '.DS_Store' in stain:
                    continue
                if refsta in stain[: -12]:
                    print('Analyzing of {} ...'.format(refsta))
                    edge_cake = ''
                    cake = ''
                    seq = ''
                    SBs = ''
                    coord_g = coord.loc[coord.genome == stain[: -12]]
                    fasta = read(args.orgs + org + '/All/' + stain)
                    for line in fasta:
                        seq += line.sequence
                    for i in coord_g.index:
                        if coord_g.start[i] < 0 and coord_g.start[i] > coord_g.end[i]:
                            seq_i = seq[coord_g.end[i]: coord_g.start[i]]
                        else:
                            seq_i = seq[coord_g.start[i]: coord_g.end[i]]
                        if args.cut[-1] == '%':
                            if len(seq_i) < 50:
                                continue
                            left_edge = seq_i[: int(len(seq_i) *  float(args.cut[: -1]) * 0.01)]
                            right_edge = seq_i[len(seq_i) -int(len(seq_i) *  float(args.cut[: -1]) * 0.01):]
                            mid = seq_i[int(len(seq_i) *  float(args.cut[: -1]) * 0.01): len(seq_i) - int(len(seq_i) *  float(args.cut[: -1]) * 0.01)]
                            cake += ('>sb{}'.format(i))
                            cake += ('\n')
                            cake += (mid)
                            cake += ('\n')
                            edge_cake += '>sb{}_right\n'.format(i)
                            edge_cake += right_edge + '\n'
                            edge_cake += '>sb{}_left\n'.format(i) 
                            edge_cake += left_edge + '\n'
                            if args.sb_out != 'NO':
                                SBs += '>sb{}\t{}\t[{}:{}]\n'.format(i, stain[: -12], coord_g.start[i], coord_g.end[i])
                                SBs += seq_i + '\n'

                        else:
                            if 0.2 * len(seq_i) < int(args.cut):
                                continue
                            left_edge = seq_i[: int(args.cut)]
                            right_edge = seq_i[len(seq_i) - int(args.cut):]
                            mid = seq_i[int(args.cut): len(seq_i) - int(args.cut)]
                            cake += ('>sb{}'.format(i))
                            cake += ('\n')
                            cake += (mid)
                            cake += ('\n')
                            edge_cake += '>sb{}_right\n'.format(i)
                            edge_cake += right_edge + '\n'
                            edge_cake += '>sb{}_left\n'.format(i) 
                            edge_cake += left_edge + '\n'
                            if args.sb_out != 'NO':
                                SBs += '>sb{}\t{}\t[{}:{}]\n'.format(i, stain[: -12], coord_g.start[i], coord_g.end[i])
                                SBs += seq_i + '\n'
                    
                    savefile = open(args.save_way + org + '/' + org + '_' + stain[: -12] + '_middle.fasta', 'w')
                    savefile.write(cake)
                    savefile.close()
                    savefile = open(args.save_way  + org + '/' + org + '_' + stain[: -12] + '_edges.fasta', 'w')
                    savefile.write(edge_cake)
                    savefile.close()
                    if args.sb_out != 'NO':
                        saveSB.write(SBs)
        saveSB.close()

    print('+ one ready org (ﾟ▽ﾟ)/')

elif args.syn_fa == 'Sibelia':
    for org in database:
        if '.DS_Store' in org:
            continue
        print('Analyzing of {} ...'.format(org))
        saveSB = open(args.sb_out + 'SB_{}.fasta'.format(org), 'w')
        #Change path if you have enother
        p = pd.read_csv(
            args.out + '/{}_results_sibelia/blocks_coords.txt'.format(org), 
            header=None, 
            names=('ind', 'size', 'description'), 
            sep='\t'
            )

        check_list = {}
        for i in range(len(p)):
            if str(p['description'][i]) != 'nan':
                if str(p['description'][i]) != 'Description':
                    check_list[str(p['description'][i])] = {}   
            else:
                break

        orgs = os.listdir(args.orgs + org + '/All/')

        for sta in orgs:
            if '.DS_Store' in sta:
                continue
            op_pa = read(args.orgs + org + '/All/{}'.format(sta))
            for line in op_pa:
                for num in check_list:
                    if num in line.description:
                        check_list[num] = [sta[: -12]]
        
        stains = []
        for ind in check_list:
            if check_list[ind][0] not in stains:
                stains.append(check_list[ind][0])

        coord = pd.read_csv(
            args.coord + org + '_results_sibelia/' + org + '.txt', 
            sep='\t',
            header=None,
            names=['sb1', 'sb2', 'genome', 'start', 'end', 'maxlen', 'contig', 'strand']
            ) 
        spike = pd.read_csv(
            args.coord + org + '_results_sibelia/' + org + '_helper.txt', 
            sep='\t',
            header=None,
            names=['genome', 'contig', 'strand']
            ) 

        checker = 0
        dataset = os.listdir(args.orgs + org + '/All/')
 
        for refsta in stains:
            for stain in dataset:
                if '.DS_Store' in stain:
                    continue
                if refsta in stain[: -12]:
                    print('Analyzing of {} ...'.format(refsta))
                    edge_cake = ''
                    cake = ''
                    SBs = ''
                    coord_g = coord.loc[coord.genome == stain[: -12]]
                    coord_spike = spike.loc[spike.genome == stain[: -12]]
                    fasta = read(args.orgs + org + '/All/' + stain)
                    for line in fasta:
                        for i in coord_g.index:
                            if line.description.split()[0] == coord_spike.contig[i]:
                                seq = line.sequence
                                if coord_spike.strand[i] == '+':
                                    seq_i = seq[coord_g.start[i]: coord_g.end[i]]
                                else:
                                    seq_i = seq[coord_g.start[i]: coord_g.end[i]]
                                    #Making conmplementary strain 3'-> 5'
                                    seq_i = seq_i.replace('A', 't')
                                    seq_i = seq_i.replace('T', 'a')
                                    seq_i = seq_i.replace('G', 'c')
                                    seq_i = seq_i.replace('C', 'g')
                                    seq_i = seq_i.upper()

                                if args.cut[-1] == '%':
                                    if len(seq_i) < 50:
                                        continue
                                    left_edge = seq_i[: int(len(seq_i) * float(args.cut[: -1]) * 0.01)]
                                    right_edge = seq_i[len(seq_i) -int(len(seq_i) * float(args.cut[: -1]) * 0.01):]
                                    mid = seq_i[int(len(seq_i) * float(args.cut[: -1]) * 0.01): len(seq_i) - int(len(seq_i) * float(args.cut[: -1]) * 0.01)]
                                    cake += ('>sb{}'.format(i))
                                    cake += ('\n')
                                    cake += (mid)
                                    cake += ('\n')
                                    edge_cake += '>sb{}_right\n'.format(i)
                                    edge_cake += right_edge + '\n'
                                    edge_cake += '>sb{}_left\n'.format(i) 
                                    edge_cake += left_edge + '\n'
                                    if args.sb_out != 'NO':
                                        SBs += '>sb{}\t{}\t{}\t[{}:{}]\n'.format(
                                            i, stain[: -12], coord_spike.contig[i], coord_g.start[i], coord_g.end[i])
                                        SBs += seq_i + '\n'
                                else:
                                    if 0.2 * len(seq_i) < int(args.cut):
                                        continue
                                    left_edge = seq_i[: int(args.cut)]
                                    right_edge = seq_i[len(seq_i) - int(args.cut):]
                                    mid = seq_i[int(args.cut): len(seq_i) - int(args.cut)]
                                    cake += ('>sb{}'.format(i))
                                    cake += ('\n')
                                    cake += (mid)
                                    cake += ('\n')
                                    edge_cake += '>sb{}_right\n'.format(i)
                                    edge_cake += right_edge + '\n'
                                    edge_cake += '>sb{}_left\n'.format(i) 
                                    edge_cake += left_edge + '\n'
                                    if args.sb_out != 'NO':
                                        SBs += '>sb{}\t{}\t{}\t[{}:{}]\n'.format(
                                            i, stain[: -12],  coord_spike.contig[i], coord_g.start[i], coord_g.end[i])
                                        SBs += seq_i + '\n'
            
                    savefile = open(args.save_way + org + '/' + org + '_' + stain[: -12] + '_middle.fasta', 'w')
                    savefile.write(cake)
                    savefile.close()
                    savefile = open(args.save_way  + org + '/' + org + '_' + stain[: -12] + '_edges.fasta', 'w')
                    savefile.write(edge_cake)
                    savefile.close()
                    if args.sb_out != 'NO':
                        saveSB.write(SBs)
        saveSB.close()

    print('+ one ready org (ﾟ▽ﾟ)/')
  
else:
    print('Error!')
                
print('Done! (´｡• ω •｡`)')