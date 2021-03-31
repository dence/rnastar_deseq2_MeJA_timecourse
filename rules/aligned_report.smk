#Daniel ence
#December 21, 2020

from scipy.stats import variation
import numpy as np

rule get_aligned_report:
	input:
		samples=np.unique(
		expand("results/star/{sample}-{unit}.Log.final.out",
				sample=units["sample"].tolist(),unit=units["unit"].tolist())).tolist()
		#samples=expand("results/star/{sample}-{unit}.Log.final.out",sample=units["sample"].tolist(),unit=units["unit"].tolist())

	output:
		"results/reports/star_percent_aligned_report.txt"
	params:
		samples=np.unique(samples["sample"].tolist()).tolist(),
		units=np.unique(units["unit"].tolist()).tolist()
	log:
		"logs/reports/star_percent_aligned_report.log"
	script:
		"../scripts/percent_aligned_report.py"

rule get_aligned_spike_report:
	input:
		expand("results/star_spike_in/{sample}-{unit}.Log.final.out",sample=units["sample"],unit=units["unit"])
	output:
		"results/reports/spike_in_percent_aligned.tsv"
	params:
		samples=samples["sample"].tolist(),
		units=units["unit"].tolist()
	log:
		"logs/reports/spike_in_star_percent_aligned_report.log"
	script:
		"../scripts/percent_aligned_report.py"

rule get_aligned_rRNA_report:
	input:
		expand("results/star_rRNA/{sample}-{unit}.Log.final.out",sample=units["sample"],unit=units["unit"])
	output:
		"results/reports/rRNA_percent_aligned_report.txt"
	params:
		samples=samples["sample"].tolist(),
		units=units["unit"].tolist()
	log:
		"logs/reports/rRNA_star_percent_aligned_report.log"
	script:
		"../scripts/percent_aligned_report.py"

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
