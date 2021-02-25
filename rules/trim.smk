
rule cutadapt_pe:
	input:
		get_fq
	output:
		fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
		fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
		qc="results/trimmed/{sample}-{unit}.qc.txt"
	params:
		"-g {} {} -G {}".format(config["trimming"]["adapter1"], config["params"]["cutadapt-pe"], config["trimming"]["adapter2"])
	log:
		"logs/cutadapt/{sample}-{unit}.log"
	shell:
		"unset TMPDIR; module load cutadapt; cutadapt {params} -o {output.fastq1} -p {output.fastq2} {input} > {output.qc}"

#rule cutadapt:
#    input:
#        get_fastq
#    output:
#        fastq="trimmed/{sample}-{unit}.fastq.gz",
#        qc="trimmed/{sample}-{unit}.qc.txt"
#    params:
#        "-a {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-se"])
#    log:
#        "logs/cutadapt/{sample}-{unit}.log"
#    wrapper:
#        "0.17.4/bio/cutadapt/se"
