
from random import shuffle
from subprocess import run
from os.path import join
import sys

'''
Sequence data for viral introns were obtained via the NCBI Virus database ('RefSeq' > 'Homo sapiens (human), taxid:9606').
The locations and sequences were extracted from the GenBank data. This dataset was filtered to remove duplicates,
resulting in a final dataset of 126 viral introns.
For each viral intron, this script uses RNAfold to fold the original sequence as well as 10 shuffled sequences
that have the same length/nucleotide composition. 
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]

	intron_file = join(project_dir, 'virus_master_with_shuffle.txt')
	out_path = join(project_dir, 'virus_introns_observed_shuffled_mfe.txt')

	introns_done, intron_set = 0, 1
	out_file = open('tmp.fa', 'w')
	with open(intron_file) as in_file, open('tmp.fa', 'w') as out_file:
		next(in_file)
		for line in in_file:
			virus, virus_score, accession_num, gene_num, _, _, _, intron_seq = line.strip().split('\t')[:8]
			intron_seq = intron_seq.upper()
			intron_info = '_'.join([virus, virus_score, accession_num, gene_num])
			intron_seq_list = list(intron_seq)
			out_file.write(f'>{intron_info}_Observed\n{intron_seq}\n')
			for i in range(1,11):
				shuffle(intron_seq_list)
				intron_seq_shuffled = ''.join(intron_seq_list)
				out_file.write(f'>{intron_info[:-3]}_Shuffled{i}\n{intron_seq_shuffled}\n')
			
	fold_call = f'RNAfold --noPS --jobs=35 --unordered --infile=tmp.fa --outfile={out_path}'
	run(fold_call.split(' '))
	run('rm tmp.fa'.split(' '))
	


