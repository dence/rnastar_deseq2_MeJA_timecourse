include: "common.smk"


rule spike_in_align:
	input:
		fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
		fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz"
	output:
		bam="results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
		star_output="results/star_spike_in/{sample}-{unit}.Log.final.out"
	log:
		"logs/spike_in_align/{sample}-{unit}.log"
	params:
		star_prefix="results/star_spike_in/{sample}-{unit}.",
		index=config["ref"]["spike_in_index"],
		extra=config["params"]["star"]
	threads:
		10
	shell:
		"unset TMPDIR; module load star/2.7.5c; set -euo pipefail;  " +
		"STAR --readFilesCommand zcat {params.extra} --runThreadN {threads} --genomeDir {params.index} " +
		"--readFilesIn {input.fastq1} {input.fastq2} --outSAMtype BAM SortedByCoordinate " +
		"--outFileNamePrefix {params.star_prefix} --outStd Log  > {log} 2>&1"

rule star_align:
	input:
		fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
		fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz"
	output:
        # see STAR manual for additional output files
		bam="results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
		star_output="results/star/{sample}-{unit}.Log.final.out"
	log:
		"logs/star/{sample}-{unit}.log"
	params:
		# path to STAR reference genome index
		star_prefix="results/star/{sample}-{unit}.",
		index=config["ref"]["index"],
        # optional parameters
		extra=config["params"]["star"]
	threads: 24
	shell:
		"unset TMPDIR; module load star/2.7.5c; set -euo pipefail;  " +
		"STAR --readFilesCommand zcat {params.extra} --runThreadN {threads} --genomeDir {params.index} " +
		"--readFilesIn {input.fastq1} {input.fastq2} --outSAMtype BAM SortedByCoordinate " +
		"--outFileNamePrefix {params.star_prefix} --outStd Log  > {log} 2>&1"
