def get_fq(wildcards):
	if config["trimming"]["skip"]:
		# no trimming, use raw reads
		return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
	else:
		# yes trimming, use trimmed data
		if not is_single_end(**wildcards):
			# paired-end sample
			return expand("results/trimmed/{sample}-{unit}.{group}.fastq.gz",
				group=[1, 2], **wildcards)
        # single end sample
		return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule make_spike_in_ref:
	input:
		ref=config["ref_fasta"],
		spike_in_ref=config["spike_in_fasta"]
	output:
		fasta="results/ref_plus_spike_in/ref_plus_spike_in.fasta"
	log:
		"logs/ref_plus_spike_in/ref_plus_spike_in.log"
	shell:
		"cat {input.ref} {input.spike_in_ref} > {output.fasta}"

rule make_ref_index:
		input:
			ref_fasta="results/ref_plus_spike_in/ref_plus_spike_in.fasta"
		output:
			"results/ref_plus_spike_in_index/Log.out"
		params:
			genome_dir="results/ref_plus_spike_in_index/"
		shell:
			"unset TMPDIR; module load star/2.7.5c; set -euo pipefail " +
			"STAR --runThreadN 1 --runMode genomeGenerate "  +
			"--genomeDir {params.genome_dir} " +
			"--genomeFastaFiles {input.ref_fasta}" +
			"--limitGenomeGenerateRAM=58220435040 --genomeSAindexNbases 11"
