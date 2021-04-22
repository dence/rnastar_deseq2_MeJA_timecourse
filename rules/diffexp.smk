def get_strandness(units):
	if "strandedness" in units.columns:
		return units["strandedness"].tolist()
	else:
		strand_list=["none"]
		return strand_list*units.shape[0]

import numpy as np

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
	shell:
		"module load htseq; htseq-count {params.option} {input.bam} {params.ref} > {output} 2> {log}"

rule htseq_count_rRNA:
	input:
		bam="results/star_rRNA/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
		bai="results/star_rRNA/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
	output:
		"results/htseq-counts_rRNA/{sample}-{unit}.htseq-counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		#ref=config["ref"]["index"]
		ref=config["ref"]["rRNA_gtf"],
		option=" --format=bam --stranded=no --type=transcript --idattr=gene_id --quiet --mode=union --nonunique=all "
	log:
		"logs/htseq-counts_rRNA/htseq-counts.{sample}-{unit}.log"
	conda:
		"../envs/htseq.yaml"
	shell:
		"module load htseq; htseq-count {params.option} {input.bam} {params.ref} > {output} 2> {log}"

rule htseq_count_spike_in:
	input:
		bam="results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
		bai="results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
	output:
		"results/htseq-counts_spike_in/{sample}-{unit}.htseq-counts.tsv"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		#ref=config["ref"]["index"]
		ref=config["ref"]["spike_in_gtf"],
		option=" --format=bam --stranded=no --type=transcript --idattr=gene_id --quiet --mode=union --nonunique=all "
	log:
		"logs/htseq-counts_rRNA/htseq-counts.{sample}-{unit}.log"
	conda:
		"../envs/htseq.yaml"
	shell:
		"module load htseq; htseq-count {params.option} {input.bam} {params.ref} > {output} 2> {log}"

#need to rewrite this rule to collate htseq-count outputs to a single count by reps matrix
rule count_matrix_with_reps:
	input:
		np.unique(expand("results/htseq-counts/{sample}-{unit}.htseq-counts.tsv", sample=units["sample"],unit=units["unit"])).tolist()
	output:
		"results/counts/counts_with_reps.tsv"
	log:
		"logs/counts/counts_with_reps.log"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		ref=config["ref"]["ref_fasta"],
		strand=get_strandness(units)
	conda:
		"../envs/htseq.yaml"
	script:
		"../scripts/count-matrix-by-reps.py"

rule spike_in_count_matrix_with_reps:
	input:
		np.unique(expand("results/htseq-counts_spike_in/{sample}-{unit}.htseq-counts.tsv", sample=units["sample"],unit=units["unit"])).tolist()
	output:
		"results/counts/spike_in_with_reps.counts.tsv"
	log:
		"logs/counts/spike_in.counts.log"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		ref=config["ref"]["spike_in_ref"],
		strand=get_strandness(units)
	conda:
		"../envs/htseq.yaml"
	script:
		"../scripts/count-matrix-by-reps.py"

rule rRNA_count_matrix_with_reps:
	input:
		np.unique(expand("results/htseq-counts_rRNA/{sample}-{unit}.htseq-counts.tsv", sample=units["sample"],unit=units["unit"])).tolist()
	output:
		"results/counts/rRNA_with_reps.counts.tsv"
	log:
		"logs/counts/rRNA.counts.log"
	params:
		samples=units["sample"].tolist(),
		units=units["unit"].tolist(),
		ref=config["ref"]["rRNA_fasta"],
		strand=get_strandness(units)
	conda:
		"../envs/htseq.yaml"
	script:
		"../scripts/count-matrix-by-reps.py"

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

rule count_matrix_rRNA:
	input:
		"results/counts/rRNA_with_reps.counts.tsv"
	output:
		"results/counts/rRNA.counts.tsv"
	params:
		sample_file=config["samples"]
	log:
		"logs/counts/rRNA_count_matrix_all.log"
	shell:
		"""
		module load R;
		Rscript scripts/count-matrix-all.R {input} {params} {output} {log}
		"""

rule count_matrix_spike_in:
	input:
		"results/counts/spike_in_with_reps.counts.tsv"
	output:
		"results/counts/spike_in.counts.tsv"
	params:
		sample_file=config["samples"]
	log:
		"logs/counts/spike_in_count_matrix_all.log"
	shell:
		"""
		module load R;
		Rscript scripts/count-matrix-all.R {input} {params} {output} {log}
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
