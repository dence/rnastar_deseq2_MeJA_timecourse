#Daniel ence
#october 6, 20202

import pysam
import pandas as pd
import numpy as np


def get_sample_counts(bam_file, reference):

	bamfile_obj = pysam.AlignmentFile(bam_file,'rb')

	ref = open(reference + "/chrName.txt",'r')

	curr_file_counts = []
	ref_list = []
	for seq in ref:
		seq = seq.strip()
		curr_seq_reads = bamfile_obj.fetch(seq)
		read_counter = 0
		for curr_read in curr_seq_reads:
			read_counter = read_counter + 1
		curr_file_counts.extend([read_counter])
		ref_list.extend([seq])

	counts = pd.DataFrame(curr_file_counts)
	counts = counts.transpose()
	counts.columns=ref_list
	counts = counts.transpose()
	return counts

############################################################


bam_files = snakemake.input.bams
counts = [get_sample_counts(f, snakemake.params.ref)
			for f in bam_files]

samples = snakemake.params.samples
print("checking sample order in the python script")
print(samples)
for t, sample in zip(counts, samples):
	t.columns = [sample]
print("checking sample order in the python script after the zip")
for t in counts:
	print(t.columns)
matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
print("Checking columns in matrix before groupby")
print(matrix.columns)
matrix = matrix.groupby(matrix.columns, axis=1).sum()
print("Checking columns in matrix")
print(matrix.columns)
matrix.to_csv(snakemake.output[0], sep="\t")
