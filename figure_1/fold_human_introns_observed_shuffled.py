
import gzip
from random import shuffle
from subprocess import run
from os.path import join
import sys

'''
The CCDS annotation of the hg19 genome was filtered to obtain human introns of length <= 5000nt
For each human intron, this script uses RNAfold to fold the original sequence as well as 10 shuffled sequences
that have the same length/nucleotide composition. 
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]

	intron_file = join(project_dir, 'hg19.gencode.basic.v45.CCDS_introns_seqs.txt.gz')
	mfe_out_prefix = 'CCDS_length_5000_introns_observed_shuffled_mfe'

	introns_done, intron_set = 0, 1
	out_file = open('tmp.fa', 'w')
	with gzip.open(intron_file, 'rt') as in_file:
		for line in in_file:
			intron_info, intron_seq = line.strip().split('\t')
			tx_id, ccds_id, gene_name, chrom, strand, intron_start, intron_end = intron_info[:-3].split('_')
			intron_start, intron_end = int(intron_start), int(intron_end)
			if intron_end-intron_start <= 5000:
				intron_seq = intron_seq.upper()
				intron_seq_list = list(intron_seq)
				out_file.write(f'>{intron_info[:-3]}_Observed\n{intron_seq}\n')
				for i in range(1,11):
					shuffle(intron_seq_list)
					intron_seq_shuffled = ''.join(intron_seq_list)
					out_file.write(f'>{intron_info[:-3]}_Shuffled{i}\n{intron_seq_shuffled}\n')
				introns_done += 1
				if introns_done % 1000 == 0:
					out_file.close()
					print(f'Done shuffling {introns_done//1000}k introns...')
					mfe_out_path = f'{mfe_out_prefix}_{intron_set}.txt'
					fold_call = f'RNAfold --noPS --jobs=35 --unordered --infile=tmp.fa --outfile={mfe_out_path}'
					run(fold_call.split(' '))
					print(f'Done folding {introns_done//1000}k introns...')
					out_file = open('tmp.fa', 'w')
					intron_set += 1
	
	out_file.close()
	print(f'Done shuffling {introns_done//1000}k introns...')
	mfe_out_path = f'{mfe_out_prefix}_{intron_set}.txt'
	fold_call = f'RNAfold --noPS --jobs=35 --unordered --infile=tmp.fa --outfile={mfe_out_path}'
	run(fold_call.split(' '))
	print(f'Done folding {introns_done//1000}k introns...')
	run('rm tmp.fa')
	


