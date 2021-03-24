rule erccdashboard:
	input:
		"results/counts/spike_in.counts.tsv"
	#output:
		#"results/spike_ins #not sure what output file type should be
	log:
		"logs/reports/erccdashboard.txt"
  	script:
  		"../scripts/spike.in.erccdashboard.R"
