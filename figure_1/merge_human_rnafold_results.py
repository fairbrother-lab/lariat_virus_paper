
import os
from numpy import mean, median
import gzip
from os.path import join
import sys

'''
This script parses the output of RNAfold for observed and shuffled human introns (see fold_human_introns_observed_shuffled.py)
and combines it into a data table for further analysis
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]
	out_path = join(project_dir, 'CCDS_length_5000_introns_observed_shuffled_mfe_merged.txt.gz')

	intron_mfe = {}
	with gzip.open(out_path, 'wt') as out_file:
		shuffled_labs = '\t'.join([f'shuffled_mfe_{i}' for i in range(1,11)])
		out_file.write('gene_name\ttranscript_id\tccds_id\tchrom\tstrand\tintron_start\tintron_end\tintron_seq\t')
		out_file.write(f'observed_mfe\t{shuffled_labs}\tmean_shuffled_mfe\tmedian_shuffled_mfe\n')
		for i in range(1,161):
			with open(join(project_dir, f'CCDS_length_5000_introns_observed_shuffled_mfe_{i}.txt')) as in_file:
				info_line, seq_line, structure_line = [in_file.readline().strip() for l in range(3)]
				while structure_line:
					tx_id, ccds_id, gene_name, chrom, strand, start, end, seq_type = info_line[1:].split('_')
					intron_id = f'{chrom}_{strand}_{start}_{end}'
					seq_mfe = structure_line.split('(')[-1][:-1]
					if intron_id not in intron_mfe:
						intron_mfe[intron_id] = {'info':[gene_name, tx_id, ccds_id, chrom, strand, start, end], 'shuffled':{}}
					if seq_type == 'Observed':
						seq_dna = seq_line.replace('U', 'T')
						intron_mfe[intron_id]['observed'] = (seq_dna, seq_mfe)
					else:
						shuffle_num = seq_type.lstrip('Shuffled')
						intron_mfe[intron_id]['shuffled'][int(shuffle_num)] = float(seq_mfe)
					
					if 'observed' in intron_mfe[intron_id] and len(intron_mfe[intron_id]['shuffled'])==10:
						intron_seq, observed_mfe = intron_mfe[intron_id]['observed']
						shuffled_mfe = [intron_mfe[intron_id]['shuffled'][i] for i in range(1,11)]
						mean_shuffled_mfe, median_shuffled_mfe = mean(shuffled_mfe), median(shuffled_mfe)
						output_str = '\t'.join([str(e) for e in intron_mfe[intron_id]['info']+[intron_seq, observed_mfe]+shuffled_mfe+[mean_shuffled_mfe, median_shuffled_mfe]])
						out_file.write(f'{output_str}\n')
						del intron_mfe[intron_id]
					
					info_line, seq_line, structure_line = [in_file.readline().strip() for l in range(3)]
			
			print(f'Done parsing MFE file {i}')

			
