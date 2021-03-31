rule erccdashboard:
	input:
		"results/counts/spike_in.counts.tsv" 
        #do we need to rbind counts.tsv? there will be an issue since this was not done by lane...
	output:
		"results/spike_ins/transgenic.TRT.CTL.pdf
        "results/spike_ins/timecourse.TRT.CTL.pdf
        #output files are 3 pdfs for each experiment
	log:
		"logs/reports/erccdashboard.txt"
  	script:
  		"../scripts/spike.in.erccdashboard.R"
