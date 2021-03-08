#Daniel Ence
#oct. 8, 2020

rule samtools_index_merged_bam:
	input:
		"results/merged_lane_bams/{sample}.merged.bam"
	output:
		"results/merged_lane_bams/{sample}.merged.bam.bai"
	log:
		"logs/merged/{sample}.samtools_index.log"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"

rule samtools_index:
	input:
		"results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam"
	output:
		"results/star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai"
	log:
		"logs/star_index.{sample}-{unit}.log"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"

rule samtools_spikein_index:
	input:
		"results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam"
	output:
		"results/star_spike_in/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai"
	log:
		"logs/star_spike_in_index.{sample}-{unit}.log"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"
