

rule RUVseq_init:
	input:
		"results/counts/counts_with_reps.tsv",
		"results/counts/spike_in.counts.tsv"
	output:
		"results/RUVseq/all_with_spike_ins.rds"
	params:
		sample_file=config["samples"]
	log:
		"logs/RUVseq/RUVseq_init.log"
	shell:
		"""
		#module load R;
		Rscript scripts/RUVseq_init.R {input[0]} {input[1]} {output} {params.sample_file} {log}
		"""

#put the rule to run RUVseq here
