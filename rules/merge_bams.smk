rule merge_lane:
    input:
        L7_bam="results/star/{sample}-lane7.Aligned.sortedByCoord.out.bam",
        L8_bam="results/star/{sample}-lane8.Aligned.sortedByCoord.out.bam"
    output:
        "results/merged_lane_bams/{sample}.merged.bam"
    log:
        "logs/picard_merge_sam_files/{sample}.log"
    shell:
        "unset TMPDIR; module load picard; java -jar $HPC_PICARD_DIR/picard.jar MergeSamFiles I={input.L7_bam} I={input.L8_bam} O={output} &> {log}"
#making a dumb assumption about the names of the bams to merged. specific to the Fr1 project. DE
#L4_bam="rmduped_reads/{sample}_{unit}.bam",
#L5_bam="rmduped_reads/{sample}_{unit}.bam"
#lambda w: expand("results/mapped/{sample}-{unit}",
#    sample = w.sample,
#    unit = units.loc[units['sample'] == w.sample].unit.to_list()
#)
