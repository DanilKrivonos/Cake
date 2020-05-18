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

parser = argparse.ArgumentParser(description='Use to find kmers')
parser.add_argument('-modifa', type=str, help='modifying_fasta', required=True)
parser.add_argument('-reffa', type=str, help='reference_fasta', required=True)
parser.add_argument('-k', type=int, help='print_kmer', default=2)
parser.add_argument('-type', type=str, help='type_of_analyze(full-1-3mers)', default='Not_full')
parser.add_argument('-save_way', type=str, help='path_of_save', required=True)
parser.add_argument('-image_save', type=str, help='image_save_path', required=True)
args = parser.parse_args()
ref_path = os.listdir(args.reffa)
k = args.k

if args.type == 'full' and k <= 3:
    while not k > 3:
        if k > 2:
            mpl.rcParams['font.size'] = 5
        if k > 1:
            mpl.rcParams['figure.figsize'] = [25, 8]
        for path in ref_path:
            if '.DS_Store' not in path:
                print('Analizing of {}'.format(path))
                print('\n')
                genome = read(args.reffa + '{}/{}_All_SB.fasta'.format(path, path))
                syntheny_bloks = read(args.modifa + '{}/{}_ENDS_SB.fasta'.format(path, path))
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
                print('\n')

                lable_k = []
                all_sb_friq = []
                end_sb_friq = []
                sort_sb_friq = []
                sort_end_sb_friq = []
                all_p_values = []
                sort_p_values = []
                sort_lable_k = []
                cake = ''
                for kmer in create_lib(k):
                    print('Counting {} ...'.format(kmer))
                    lable_k.append(kmer)
                    all_sb_friq.append(new_genome.count(kmer) * 100/(len(new_genome) - k + 1))
                    end_sb_friq.append(all_ends.count(kmer) * 100/(len(all_ends) - k + 1))
                    cake += 'Frequency in ends {} \t{}\n'.format(kmer, str(all_ends.count(kmer) * 100/(len(all_ends) - k + 1)))
                    cake += 'Frequency in genome {} \t{}\n'.format(kmer, str(new_genome.count(kmer) * 100/(len(new_genome) - k + 1)))
                    cake += 'Ratio \t{}\n'.format(str((all_ends.count(kmer) * 100/(len(all_ends) - k + 1))/((new_genome.count(kmer) * 100)/(len(new_genome) - k + 1))))
                    pvalue = stats.chi2_contingency([
                        [len(all_ends), len(new_genome)],
                        [all_ends.count(kmer), new_genome.count(kmer)]
                        ])[1]
                    print('P-value = {}'.format(pvalue))
                    print('\n')
                    cake += 'P-value for {} \t{}\n'.format(kmer, str(math.log1p(pvalue)))
                    all_p_values.append(pvalue)
                
                #correction for multiple comparisons
                all_p_values = multipletests(all_p_values, method='fdr_bh')[1]
                all_p_values = list(np.array(all_p_values))
                sort_p_values = list(np.sort(np.array(all_p_values)))
                
                for sornum in sort_p_values:
                    for num in range(len(all_p_values)):
                        if sornum == all_p_values[num]:
                            #0 добавляется дважды, что если пи валбю равны
                            if all_sb_friq[num] not in sort_sb_friq:
                                sort_sb_friq.append(all_sb_friq[num])
                                sort_end_sb_friq.append(end_sb_friq[num])
                                sort_lable_k.append(lable_k[num])
                          
                                     
                savefile = open(args.save_way + '{}/{}_{}mers.txt'.format(path, path, k), 'w')
                savefile.write(cake)
                savefile.close()
                
                print('Text file is complete!\n༼ つ ◕_◕ ༽つ')
                print('\n')
                print('Now please wait a picture\n(ﾉ◕ヮ◕)ﾉ*:･ﾟ✧')
                print('\n')

                x = np.arange(len(sort_lable_k))
                width = 0.35
                fig, ax = plt.subplots()
                rects1 = ax.bar(x - width/2, sort_sb_friq, width, label='K-mer in all SB, %')
                rects2 = ax.bar(x + width/2, sort_end_sb_friq, width, label='K-mer Ends, %')

                ax.set_ylabel('Frequency, %')
                ax.set_xlabel('K-mers')
                ax.set_title('{}'.format(path))
                ax.set_xticks(x)
                if k < 3:
                    ax.set_xticklabels(sort_lable_k)
                else:
                    ax.set_xticklabels(sort_lable_k, rotation='vertical')
        
                if k < 2:
                    ax.legend(loc='center', bbox_to_anchor=(1.15, 0.5))
                else:
                    ax.legend(loc='center', bbox_to_anchor=(1.05, 0.5))

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
                plt.savefig(args.image_save + '{}/{}_{}mers.pdf'.format(path, path, k), format='pdf', dpi=300)
                #plt.show()

                print('\n')
                print('Analyze of {} was done!'.format(path))
                print('\n')

        print('{}-mers analyze was done'.format(k))
        print('\n')
        k += 1

else:
    if k > 1:
        mpl.rcParams['figure.figsize'] = [25, 8]
    if k > 2:
        mpl.rcParams['font.size'] = 5
    for path in ref_path:
        if '.DS_Store' not in path:
            print('Analizing of {}'.format(path))
            print('\n')
            genome = read(args.reffa + '{}/{}_All_SB.fasta'.format(path, path))
            syntheny_bloks = read(args.modifa + '{}/{}_ENDS_SB.fasta'.format(path, path))
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
            print('\n')

            lable_k = []
            all_sb_friq = []
            end_sb_friq = []
            sort_sb_friq = []
            sort_end_sb_friq = []
            all_p_values = []
            sort_p_values = []
            cake = ''
            for kmer in create_lib(k):
                print('Counting {} ...'.format(kmer))
                lable_k.append(kmer)
                all_sb_friq.append(new_genome.count(kmer) * 100/(len(new_genome) - k + 1))
                end_sb_friq.append(all_ends.count(kmer) * 100/(len(all_ends) - k + 1))
                cake += 'Frequency in ends {} \t{}\n'.format(kmer, str(all_ends.count(kmer) * 100/(len(all_ends) - k + 1)))
                cake += 'Frequency in genome {} \t{}\n'.format(kmer, str(new_genome.count(kmer) * 100/(len(new_genome) - k + 1)))
                cake += 'Ratio \t{}\n'.format(str((all_ends.count(kmer) * 100/(len(all_ends) - k + 1))/((new_genome.count(kmer) * 100)/(len(new_genome) - k + 1))))
                pvalue = stats.chi2_contingency([
                    [len(all_ends), len(new_genome)],
                    [all_ends.count(kmer), new_genome.count(kmer)]
                    ])[1]
                print('P-value = {}'.format(pvalue))
                print('\n')
                cake += 'P-value for {} \t{}\n'.format(kmer, str(math.log1p(pvalue)))
                all_p_values.append(pvalue)

            #correction for multiple comparisons
            all_p_values = multipletests(all_p_values, method='fdr_bh')[1]
            all_p_values = list(np.array(all_p_values))
            sort_p_values = list(np.sort(np.array(all_p_values)))

            for sornum in sort_p_values:
                    for num in range(len(all_p_values)):
                        if sornum == all_p_values[num]:
                            #0 добавляется дважды, что если пи валбю равны
                            if all_sb_friq[num] not in sort_sb_friq:
                                sort_sb_friq.append(all_sb_friq[num])
                                sort_end_sb_friq.append(end_sb_friq[num])
                                sort_lable_k.append(lable_k[num])           
            savefile = open(args.save_way + '{}/{}_{}mers.txt'.format(path, path, k), 'w')
            savefile.write(cake)
            savefile.close()

            print('Text file is complete!\n༼ つ ◕_◕ ༽つ')
            print('\n')
            print('Now please wait a picture\n(ﾉ◕ヮ◕)ﾉ*:･ﾟ✧')
            print('\n')

            x = np.arange(len(lable_k))
            width = 0.35
            fig, ax = plt.subplots()
            rects1 = ax.bar(x - width/2, sort_sb_friq, width, label='K-mer in all SB, %')
            rects2 = ax.bar(x + width/2, sort_end_sb_friq, width, label='K-mer Ends, %')

            ax.set_ylabel('Frequency, %')
            ax.set_xlabel('K-mers')
            ax.set_title('{}'.format(path))
            ax.set_xticks(x)
            if k < 3:
                ax.set_xticklabels(lable_k)
            else:
                ax.set_xticklabels(lable_k, rotation='vertical', fontsize='7')
            if k < 2:
                ax.legend(loc='center', bbox_to_anchor=(1, 0.5))
            else:
                ax.legend(loc='center', bbox_to_anchor=(1.05, 0.5))

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
            autolabel(rects2)
            fig.tight_layout()
            plt.savefig(args.image_save + '{}/{}_{}mers.pdf'.format(path, path, k), format='pdf', dpi=300, )
            #plt.show()

            print('\n')
            print('Analyze of {} was done!'.format(path))
            print('\n')

    print('{}-mers analyze was done'.format(k))
    print('\n')