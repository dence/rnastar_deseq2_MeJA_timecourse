#Daniel ence
#october 6, 20202

import pysam
import pandas as pd
import numpy as np
from scipy.stats import variation
import re

def get_number_of_reads_in_cat(star_log_file, target_string):
	star_log = open(star_log_file,'r')
	for line in star_log:
		if(target_string in line):
			regex = r"\s\|\s(\d+)"
			match_list = re.findall(regex, line, flags=0)
			return match_list[0]

def get_sample_counts(star_log_file):

	#print("in get sample counts")
	#print("this is star log file")
	#print(star_log_file)

	curr_file_counts = []
	target_strings = ["Number of input reads",
	"Uniquely mapped reads number",
	"Number of reads mapped to multiple loci",
	"Number of reads mapped to too many loci",
	"Number of reads unmapped: too many mismatches",
	"Number of reads unmapped: too short",
	"Number of reads unmapped: other"]

	for string in target_strings:
		curr_file_counts.append(
			get_number_of_reads_in_cat(star_log_file,string))

	counts = pd.DataFrame(curr_file_counts)
	counts = counts.transpose()
	counts.columns=target_strings
	counts = counts.transpose()

	return counts

############################################################
star_logs = np.unique(snakemake.input)
#print("Passing this to get_sample_counts")
#print(star_logs)
#print(len(star_logs))
counts = [get_sample_counts(f)
			for f in star_logs]

samples = snakemake.params.samples
#print("checking sample order in the python script")
#print(samples)
for t, sample in zip(counts, samples):
	t.columns = [sample]
#print("checking sample order in the python script after the zip")
#for t in counts:
#	print(t.columns)
matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
#print(matrix.columns)
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix = matrix.transpose()
#print("Checking columns in matrix")
#print(matrix.columns)
matrix["proportion_of_library_reads_out_of_all_reads"] = \
	matrix["Number of input reads"] / matrix["Number of input reads"].sum()
matrix["percent_of_library_reads_mapped"] = \
	(matrix["Uniquely mapped reads number"] + matrix["Number of reads mapped to multiple loci"]) / matrix["Number of input reads"]
	
matrix.to_csv(snakemake.output[0], sep="\t")
