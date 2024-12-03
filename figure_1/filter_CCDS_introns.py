
import gzip
from os.path import join
import sys

'''
This script takes the GENCODE basic annotation and outputs a BED file
containing introns from protein-coding, CCDS-annotated genes
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]

	anno_file = join(project_dir, 'hg19.gencode.basic.v45.annotation.gff3.gz')
	out_path = join(project_dir, 'hg19.gencode.basic.v45.CCDS_introns.bed.gz')

	gene_tx = {}
	with gzip.open(anno_file, 'rt') as in_file:
		for line in in_file:
			if line[0] != '#':
				chrom, source, feature, start, end, _, strand, _, info = line.strip().split('\t')
				info = dict([tuple(e.split('=')) for e in info.split(';')])
				if info['gene_type'] == 'protein_coding':
					if 'tag' in info and 'CCDS' in info['tag'].split(','):
						if feature == 'transcript':
							gene_tx[info['transcript_id']] = {'gene_name':info['gene_name'], 'ccds_id':info['ccdsid'], 'coords':(chrom, strand, int(start), int(end)), 'exon':[], 'fivep_UTR':[], 'threep_UTR':[]}
						if feature == 'CDS':
							gene_tx[info['transcript_id']]['exon'].append((int(start), int(end)))
						if feature == 'five_prime_UTR':
							gene_tx[info['transcript_id']]['fivep_UTR'].append((int(start), int(end)))
						if feature == 'three_prime_UTR':
							gene_tx[info['transcript_id']]['threep_UTR'].append((int(start), int(end)))
	print(f'Done parsing {len(gene_tx)} transcripts...')

	tx_done = 0
	introns_done = set()
	with gzip.open(out_path, 'wt') as out_file:
		for tx_id in gene_tx:
			ccds_id = gene_tx[tx_id]['ccds_id']
			gene_name = gene_tx[tx_id]['gene_name']
			chrom, strand, tx_start, tx_end = gene_tx[tx_id]['coords']
			if strand == '-':
				gene_tx[tx_id]['exon'].reverse()
				gene_tx[tx_id]['fivep_UTR'].reverse()
				gene_tx[tx_id]['threep_UTR'].reverse()
			for gene_region in ['exon', 'fivep_UTR', 'threep_UTR']:
				for i in range(len(gene_tx[tx_id][gene_region])):
					if i < len(gene_tx[tx_id][gene_region])-1:
						intron_start, intron_end = gene_tx[tx_id][gene_region][i][1], gene_tx[tx_id][gene_region][i+1][0]-1
						intron_id = f'{chrom}_{strand}_{intron_start}_{intron_end}'
						if intron_id not in introns_done:
							intron_info = f'{tx_id.split("_")[0]}_{ccds_id}_{gene_name}_{intron_id}'
							out_file.write(f'{chrom}\t{intron_start}\t{intron_end}\t{intron_info}\t0\t{strand}\n')
							introns_done.add(intron_id)
	
			tx_done += 1
			if tx_done % 5000 == 0:
				print(f'Done processing {tx_done//1000}k transcripts...')
	