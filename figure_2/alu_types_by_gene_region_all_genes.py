
from intervaltree import Interval, IntervalTree
import gzip
from collections import Counter
from os.path import join
import sys

'''
This script takes the Repeatmasker and GENCODE basic annotations and outputs
a data table containing information on the number and type of Alu insertions 
in each gene region (3' UTR, exon/CDS, intron and 5' UTR)
'''

if __name__ == '__main__':

	project_dir = sys.argv[1]
	repeatmasker_file_hg19 = join(project_dir, 'hg19.repeatmasker.bed.gz')
	anno_file_hg19 = join(project_dir, 'hg19.gencode.basic.v45.annotation.gff3.gz')
	table_out_path = join(project_dir, 'hg19.gencode.basic.v45.gene_regions_with_alu_content.txt.gz')
	intron_bed_out_path = join(project_dir, 'hg38.gencode.basic.v45.introns_with_alu_category.bed.gz')
	
	alus, total_alu_count = {}, 0
	with gzip.open(repeatmasker_file_hg19, 'rt') as in_file:
		for line in in_file:
			chrom, start, end, info, _, strand = line.strip().split('\t')
			if 'Alu' in info:
				if chrom not in alus:
					alus[chrom] = IntervalTree()
				alus[chrom].add(Interval(int(start)+1, int(end), {'strand':strand}))
				total_alu_count += 1
	print('Done parsing Alus...')
	
	gene_tx, coding_tx = {}, set()
	with gzip.open(anno_file_hg19, 'rt') as in_file:
		for line in in_file:
			if line[0] != '#':
				chrom, source, feature, start, end, _, strand, _, info = line.strip().split('\t')
				info = dict([tuple(e.split('=')) for e in info.split(';')])
				if feature == 'transcript':
					gene_tx[info['transcript_id']] = {'info':[info['gene_name'], info['gene_id'], info['gene_type']], 'coords':(chrom, strand, int(start), int(end)), 'exon':[], 'CDS':[], 'fivep_UTR':[], 'threep_UTR':[]}
				if feature == 'CDS':
					gene_tx[info['transcript_id']]['CDS'].append((int(start), int(end)))
					coding_tx.add(info['transcript_id'])
				if feature == 'exon':
					gene_tx[info['transcript_id']]['exon'].append((int(start), int(end)))
				if feature == 'five_prime_UTR':
					gene_tx[info['transcript_id']]['fivep_UTR'].append((int(start), int(end)))
				if feature == 'three_prime_UTR':
					gene_tx[info['transcript_id']]['threep_UTR'].append((int(start), int(end)))
	print(f'Done parsing {len(gene_tx)} transcripts...')
	
	tx_done = 0
	intron_alu_set = set()
	with gzip.open(table_out_path, 'wt') as table_out, gzip.open(intron_bed_out_path, 'wt') as bed_out:
		table_out.write('gene_name\tgene_id\tgene_type\ttranscript_id\tgene_region\tchrom\tstrand\tstart\tend\ttotal_alu_count\ttotal_alu_length\tir_alu_pair_count\thas_alu\thas_ir_alu\talu_category\n')
		for tx_id in gene_tx:
			chrom, strand, tx_start, tx_end = gene_tx[tx_id]['coords']
			gene_name, gene_id, gene_type = gene_tx[tx_id]['info']
			if chrom not in alus:
				continue
			if strand == '-':
				gene_tx[tx_id]['exon'].reverse()
				gene_tx[tx_id]['CDS'].reverse()
				gene_tx[tx_id]['fivep_UTR'].reverse()
				gene_tx[tx_id]['threep_UTR'].reverse()
			gene_regions = ['CDS', 'fivep_UTR', 'threep_UTR'] if tx_id in coding_tx else ['exon', 'fivep_UTR', 'threep_UTR']
			for gene_region in gene_regions:
				for i in range(len(gene_tx[tx_id][gene_region])):
					region_start, region_end = gene_tx[tx_id][gene_region][i]
					region_alus = [a for a in alus[chrom].overlap(region_start, region_end) if a.begin>=region_start and a.end<=region_end]
					region_alu_length = len(set(sum((list(range(a.begin, a.end+1)) for a in region_alus), start=[])))
					region_total_alus = len(region_alus)
					region_alu_orientations = Counter([a.data['strand'] for a in region_alus])
					region_ir_alu_pairs = min(region_alu_orientations['+'], region_alu_orientations['-']) if len(region_alu_orientations) > 1 else 0
					has_alu, has_ir_alu = region_total_alus > 0, region_ir_alu_pairs > 0
					if has_ir_alu:
						alu_category = 'ir_alu'
					elif region_total_alus > 1:
						alu_category = 'non_ir_multi_alu'
					elif region_total_alus == 1:
						alu_category = 'non_ir_single_alu'
					else:
						alu_category = 'no_alu'
					table_out.write(f'{gene_name}\t{gene_id}\t{gene_type}\t{tx_id}\t{gene_region}\t{chrom}\t{strand}\t{region_start}\t{region_end}\t{region_total_alus}\t{region_alu_length}\t{region_ir_alu_pairs}\t{has_alu}\t{has_ir_alu}\t{alu_category}\n')
					if i < len(gene_tx[tx_id][gene_region])-1:
						intron_start, intron_end = gene_tx[tx_id][gene_region][i][1]+1, gene_tx[tx_id][gene_region][i+1][0]-1
						intron_alus = [a for a in alus[chrom].overlap(intron_start, intron_end) if a.begin>=intron_start and a.end<=intron_end]
						for alu_int in intron_alus:
							intron_alu_set.add(f'{chrom}_{alu_int.begin}_{alu_int.end}')
						intron_alu_length = len(set(sum((list(range(a.begin, a.end+1)) for a in intron_alus), start=[])))
						intron_total_alus = len(intron_alus)
						intron_alu_orientations = Counter([a.data['strand'] for a in intron_alus])
						intron_ir_alu_pairs = min(intron_alu_orientations['+'], intron_alu_orientations['-']) if len(intron_alu_orientations) > 1 else 0
						has_alu, has_ir_alu = intron_total_alus > 0, intron_ir_alu_pairs > 0
						if has_ir_alu:
							alu_category = 'ir_alu'
						elif intron_total_alus > 1:
							alu_category = 'non_ir_multi_alu'
						elif intron_total_alus == 1:
							alu_category = 'non_ir_single_alu'
						else:
							alu_category = 'no_alu'
						table_out.write(f'{gene_name}\t{gene_id}\t{gene_type}\t{tx_id}\tintron\t{chrom}\t{strand}\t{intron_start}\t{intron_end}\t{intron_total_alus}\t{intron_alu_length}\t{intron_ir_alu_pairs}\t{has_alu}\t{has_ir_alu}\t{alu_category}\n')
						bed_out.write(f'{chrom}\t{intron_start-1}\t{intron_end}\t{tx_id}-{alu_category}\t{intron_total_alus}\t{strand}\n')

			tx_done += 1
			if tx_done % 5000 == 0:
				print(f'Done processing {tx_done//1000}k transcripts...')
	
	print(f'Of {total_alu_count} annotated Alu, {len(intron_alu_set)} ({round(100*float(len(intron_alu_set))/total_alu_count, 2)}%) occur in introns')
			





				


	

				