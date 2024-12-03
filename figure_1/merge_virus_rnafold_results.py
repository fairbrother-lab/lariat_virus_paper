
from numpy import mean, median
from os.path import join
import sys

'''
This script parses the output of RNAfold for observed and shuffled viral introns (see fold_virus_introns_observed_shuffled.py)
and combines it into a data table for further analysis
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]

	fold_output = join(project_dir, 'virus_introns_observed_shuffled_mfe.txt')
	out_path = join(project_dir, 'virus_introns_observed_shuffled_mfe_merged.txt')

	intron_mfe = {}
	with open(out_path, 'w') as out_file, open(fold_output) as in_file:
		shuffled_labs = '\t'.join([f'shuffled_mfe_{i}' for i in range(1,11)])
		out_file.write('virus\tvirus_score\taccession_num\tgene_num\tintron_seq\t')
		out_file.write(f'observed_mfe\t{shuffled_labs}\tmean_shuffled_mfe\tmedian_shuffled_mfe\n')
		info_line, seq_line, structure_line = [in_file.readline().strip() for l in range(3)]
		while structure_line:
			virus, virus_score, accession_str, accession_num, gene_num, seq_type = info_line[1:].split('_')
			seq_mfe = structure_line.split('(')[-1][:-1]
			intron_id = (virus, virus_score, f'{accession_str}_{accession_num}')
			if intron_id not in intron_mfe:
				intron_mfe[intron_id] = {'shuffled':{}}
			if seq_type == 'Observed':
				seq_dna = seq_line.replace('U', 'T')
				intron_mfe[intron_id]['observed'] = (seq_dna, seq_mfe)
				intron_mfe[intron_id]['info'] = (virus, virus_score, f'{accession_str}_{accession_num}', gene_num)
			else:
				shuffle_num = seq_type.lstrip('Shuffled')
				intron_mfe[intron_id]['shuffled'][int(shuffle_num)] = float(seq_mfe)
			
			if 'observed' in intron_mfe[intron_id] and len(intron_mfe[intron_id]['shuffled'])==10:
				intron_seq, observed_mfe = intron_mfe[intron_id]['observed']
				shuffled_mfe = [intron_mfe[intron_id]['shuffled'][i] for i in range(1,11)]
				mean_shuffled_mfe, median_shuffled_mfe = mean(shuffled_mfe), median(shuffled_mfe)
				output_str = '\t'.join([str(e) for e in list(intron_mfe[intron_id]['info'])+[intron_seq, observed_mfe]+shuffled_mfe+[mean_shuffled_mfe, median_shuffled_mfe]])
				out_file.write(f'{output_str}\n')
			
			info_line, seq_line, structure_line = [in_file.readline().strip() for l in range(3)]
	
			
