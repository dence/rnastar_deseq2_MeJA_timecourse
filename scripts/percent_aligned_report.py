#Daniel ence
#october 6, 20202

import pysam
import pandas as pd
import numpy as np
from scipy.stats import variation
import re
import logging

def get_number_of_reads_in_cat(star_log_file, target_string):
	star_log = open(star_log_file,'r')
	for line in star_log:
		if(target_string in line):
			regex = r"\s\|\s(\d+)"
			match_list = re.findall(regex, line, flags=0)
			return match_list[0]

def get_sample_counts(star_log_file, log_file_handle):

	log_file_handle.write("in get sample counts\n")
	log_file_handle.write("this is star log file\n")
	log_file_handle.write(star_log_file)
	log_file_handle.write("\n")

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
	print("This is counts at the end of get_sample_counts")
	print(str(counts))
	print("\n")
	return counts

############################################################

star_logs = snakemake.input
script_log = snakemake.log
#print("What is log?\t" + str(script_log))
log_file = open(str(script_log),"w")

log_file.write("Passing this to get_sample_counts")
log_file.write(str(star_logs))
log_file.write("\n")
log_file.write(str(len(star_logs)))
log_file.write("\n")
counts = [get_sample_counts(f, log_file)
			for f in star_logs]

samples = snakemake.params.samples
units = snakemake.params.units
log_file.write("checking sample order in the python script")
log_file.write(str(samples))
for t, sample in zip(counts, star_logs):
	t.columns = [sample]
log_file.write("\nchecking sample order in the python script after the zip")
for t in counts:
	log_file.write(str(t.columns))
matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
log_file.write(str(matrix.columns))
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix = matrix.transpose()
log_file.write("Checking columns in matrix")
log_file.write(str(matrix.columns))
matrix["proportion_of_library_reads_out_of_all_reads"] = \
	matrix["Number of input reads"] / matrix["Number of input reads"].sum()
matrix["proportion_of_library_reads_mapped"] = \
	(matrix["Uniquely mapped reads number"] + matrix["Number of reads mapped to multiple loci"]) / matrix["Number of input reads"]
matrix["sensitivity (uniquely mapped / (uniquely mapped + too many mismatches))"] = \
	(matrix["Uniquely mapped reads number"] / (matrix["Uniquely mapped reads number"] + matrix["Number of reads unmapped: too many mismatches"]))
matrix["specificity (unmapped:other / (ummapped:other + too many loci))"] = \
	(matrix["Number of reads unmapped: other"] / (matrix["Number of reads unmapped: other"] + matrix["Number of reads mapped to too many loci"]))
matrix["ratio mapped unique to multiple mapped"] = \
	(matrix["Uniquely mapped reads number"] / matrix["Number of reads mapped to multiple loci"])

matrix.to_csv(snakemake.output[0], sep="\t")
