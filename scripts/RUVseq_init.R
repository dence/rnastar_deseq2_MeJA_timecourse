# Make a Merged Counts Table

# List of inputs you need:

# spike_in.counts <- # need input table with ERCC info by lane
# counts_with_reps <- # need input table with gene info by lane
# list.of.sample.names <- # need input list of sample names as character vector with length 135

args=commandArgs(trailingOnly = TRUE)
counts_with_reps_file = args[1] #tab-delim with gene counts per sample per lane
spike_in_with_reps_file = args[2] #tab-delim with spike-in counts per sample per lane
out_filename = args[3] #name of file to write final master counts table too. 
samples_table = args[4] #samples.tsv file to get a list of sample names
log_file = args[5]

log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

#Will hard code files for now, but will be reading from argv in the final
counts_with_reps <- read.table(counts_with_reps_file,header=T,sep="\t")
spike_in.counts <- read.table(spike_in_with_reps_file,header=T,sep="\t")

# First, rbind the tables together
print("counts_with_reps")
print(dim(counts_with_reps))
print("spike_in.counts")
print(spike_in.counts)
all.counts <- rbind(spike_in.counts,counts_with_reps)

# Now sum every three rows
master <- as.data.frame(sapply(seq(2,dim(all.counts)[2],by=3),function(i) rowSums(all.counts[,i:(i+2)])))

# We removed the first column, "genes", from the last function, and we want it as our rownames
rownames(master) <- all.counts$gene

saveRDS(master,file=out_filename)

# Label "master" with list of sample names \

# names(master) <- list.of.sample.names

# Split the experiments 
#transgenic.exp <- master[,1:47]
#timecourse.exp <- master[,48:135]

# Write output files for each exp, and you're done
