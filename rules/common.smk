import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
#min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

#print("here")
#print("|".join(samples.index))
#print(units.head())
#print("|".join(units["unit"]))
#print(units["unit"])

wildcard_constraints:
	sample="|".join(samples.index),
	unit="|".join(units["unit"])

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

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
