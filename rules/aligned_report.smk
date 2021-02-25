#Daniel ence
#December 21, 2020

from scipy.stats import variation

rule get_aligned_report:
	input:
		expand("results/star/{sample}-{unit}.Log.final.out",sample=units["sample"],unit=units["unit"])
	output:
		"results/reports/star_percent_aligned_report.txt"
	params:
		samples=samples["sample"].tolist()
	script:
		"../scripts/percent_aligned_report.py"

rule get_aligned_spike_report:
	input:
		expand("results/star_spike_in/{sample}-{unit}.Log.final.out",sample=units["sample"],unit=units["unit"])
	output:
		"results/reports/spike_in_aligned_report.txt"
	params:
		samples=samples["sample"].tolist()
	script:
		"../scripts/percent_aligned_report.py"

rule count_matrix_spike_in:
	input:
		bams=expand("results/star_spike_in/{unit.sample}-{unit.unit}.Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
		bai=expand("results/star_spike_in/{unit.sample}-{unit.unit}.Aligned.sortedByCoord.out.bam.bai", unit=units.itertuples())
	output:
		"results/counts/spike_in.tsv"
	params:
		samples=units["sample"].tolist(),
		ref=config["ref"]["spike_in_index"]
	conda:
		"../envs/pandas.yaml"
	script:
		"../scripts/count-matrix-bams.py"

rule get_spike_in_coefficient_of_deviation_report:
	input:
		"results/counts/spike_in.tsv"
	output:
		"results/reports/spike_in_coefficient_of_deviation.txt"
	run:
		print("Input is:\t")
		print(input[0])
		spike_in_counts = pd.read_table(input[0]).set_index("gene",drop=True)
		#foreach sample,
		#get the CoD for all the genes
		sample_names = spike_in_counts.columns
		coeff_of_deviations = []
		for curr_sample in sample_names:
			coeff_of_deviations.append(variation(spike_in_counts[curr_sample]))
		df = pd.DataFrame({ "sample": sample_names,
							"coeff_of_deviations": coeff_of_deviations })
		df.to_csv(output[0],sep="\t")
