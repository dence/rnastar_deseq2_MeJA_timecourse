def get_strandness(units):
	if "strandedness" in units.columns:
		return units["strandedness"].tolist()
	else:
		strand_list=["none"]
		return strand_list*units.shape[0]

import numpy as np

#rule count_matrix:
#	input:
#		bams=expand("results/merged_lane_bams/{sample}.merged.bam", sample=np.unique(units["sample"]).tolist()),
#		bai=expand("results/merged_lane_bams/{sample}.merged.bam.bai", sample=np.unique(units["sample"]).tolist())
#	output:
#		"results/counts/all.tsv"
#	params:
#		samples=units["sample"].tolist(),
#		ref=config["ref"]["index"]
#	log:
#		"logs/counts/count_matrix.log"
#	conda:
#		"../envs/pandas.yaml"
#	script:
#		"../scripts/count-matrix-bams.py"

rule count_matrix_spike_in:
	input:
		bams=np.unique(
		expand("results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
				sample=units["sample"].tolist(),unit=units["unit"].tolist())).tolist(),
		#bams=expand("results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam", sample=units["sample"],unit=units["unit"]),
		bai=np.unique(
		expand("results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
				sample=units["sample"].tolist(),unit=units["unit"].tolist())).tolist()
		#bai=expand("results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai", sample=units["sample"],unit=units["unit"])
	output:
		"results/counts/spike_in.counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		ref=config["ref"]["spike_in_index"]
	log:
		"logs/counts/spike_in.count_matrix.log"
	conda:
		"../envs/pandas.yaml"
	script:
		"../scripts/count-matrix-bams.py"

rule count_matrix_rRNA:
	input:
		bams=np.unique(expand("results/star_rRNA/{sample}-{unit}.Aligned.sortedByCoord.out.bam", sample=units["sample"],unit=units["unit"])).tolist(),
		bai=np.unique(expand("results/star_rRNA/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai", sample=units["sample"],unit=units["unit"])).tolist()
	output:
		"results/counts/rRNA.counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		ref=config["ref"]["rRNA_index"]
	log:
		"logs/counts/rRNA.count_matrix.log"
	conda:
		"../envs/pandas.yaml"
	script:
		"../scripts/count-matrix-bams.py"

rule htseq_count:
	input:
		bam="results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
		bai="results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
	output:
		"results/htseq-counts/{sample}-{unit}.htseq-counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		#ref=config["ref"]["index"]
		ref=config["ref"]["ref_gtf"],
		option=" --format=bam --stranded=no --type=transcript --idattr=gene_id --quiet --mode=union --nonunique=all "
	log:
		"logs/htseq-counts/htseq-counts.{sample}-{unit}.log"
	conda:
		"../envs/htseq.yaml"
	#script:
	#	"../scripts/count-matrix-bams.py"
	shell:
		"""
		module load htseq;
		htseq-count {params.option} {input.bam} {params.ref} > {output} 2> {log}
		"""

#need to rewrite this rule to collate htseq-count outputs to a single count by reps matrix
rule count_matrix_with_reps:
	input:
		bams=np.unique(expand("results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam", sample=units["sample"],unit=units["unit"])).tolist(),
		bai=np.unique(expand("results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai", sample=units["sample"],unit=units["unit"])).tolist()
	output:
		"results/htseq-counts/{sample}-{unit}.htseq-counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		#ref=config["ref"]["index"]
		ref=config["ref"]["ref_gtf"],
		option=" --stranded=no --type=transcript --quiet --mode=union --nonunique=all "
	threads: 10
	log:
		"logs/htseq-counts/htseq-counts.{sample}-{unit}.log"
	conda:
		"../envs/htseq.yaml"
	#script:
		"../scripts/count-matrix-bams.py"
	shell:
		"""
		module load htseq;
		htseq-count --nprocesses={threads} {params.option} {input} {params.ref} > {output} 2> {log}
		"""

def get_deseq2_threads(wildcards=None):
	# https://twitter.com/mikelove/status/918770188568363008
	few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
	return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
	input:
		counts="results/counts/all.tsv"
	output:
		"results/deseq2/all.rds"
	params:
		samples=config["samples"]
	conda:
		"../envs/deseq2.yaml"
	log:
		"logs/deseq2/init.log"
	threads: get_deseq2_threads()
	#script:
	#	"../scripts/deseq2-init.R"
	shell:
		"""
		module load R;
		Rscript ../scripts/deseq2-init.R {input.counts} {output} {params.samples} {log} {threads}
		"""

rule count_matrix:
	input:
		"results/counts/counts_with_reps.tsv"
	output:
		"results/counts/all.tsv"
	params:
		sample_file=config["samples"]
	log:
		"logs/counts/count_matrix_all.log"
	shell:
		"""
		module load R;
		Rscript scripts/count-matrix-all.R {input} {params} {output} {log}
		"""

rule pca:
	input:
		"results/deseq2/all.rds"
	output:
		report("results/pca.svg", "../report/pca.rst")
	params:
		pca_labels=config["pca"]["labels"]
	conda:
		"../envs/deseq2.yaml"
	log:
		"logs/pca.log"
	script:
		"../scripts/plot-pca.R"

def get_contrast(wildcards):
	return config["diffexp"]["contrasts"][wildcards.contrast]

rule deseq2:
	input:
		"results/deseq2/all.rds"
	output:
		table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
		ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
	params:
		contrast=get_contrast
	conda:
		"../envs/deseq2.yaml"
	log:
		"logs/deseq2/{contrast}.diffexp.log"
	threads: get_deseq2_threads
	script:
		"../scripts/deseq2.R"
