import pandas as pd

def get_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 1 #non stranded protocol
    elif strandedness == "yes":
        return 2 #3rd column
    elif strandedness == "reverse":
        return 3 #4th column, usually for Illumina truseq
    else:
        raise ValueError(("'strandedness' column should be empty or have the "
                          "value 'none', 'yes' or 'reverse', instead has the "
                          "value {}").format(repr(strandedness)))

script_log = snakemake.log
log_file = open(str(script_log),"w")

counts = [pd.read_table(f, index_col=0, usecols=[0, get_column(strandedness)],
          header=None, skiprows=4)
          for f, strandedness in zip(snakemake.input, snakemake.params.strand)]

for t, sample, unit in zip(counts, snakemake.params.samples, snakemake.params.units):
	t.columns = [sample + "-" + unit]
	log_file.write("Made this sample-unit pair")
	log_file.write(sample + "_" + unit)

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
