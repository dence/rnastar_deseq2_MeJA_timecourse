# Make a Merged Counts Table

# List of inputs you need:

# spike_in.counts <- # need input table with ERCC info by lane
# counts_with_reps <- # need input table with gene info by lane
# list.of.sample.names <- # need input list of sample names as character vector with length 135

args=commandArgs(trailingOnly = TRUE)
counts_with_reps_file = args[1] #tab-delim with gene counts per sample per lane
samples_file = args[2] #samples.tsv file to get a list of sample names
out_filename = args[3] #name of file to write final master counts table too. 
log_file = args[4]


#Will hard code files for now, but will be reading from argv in the final
counts_with_reps <- read.table(counts_with_reps_file,header=T,sep="\t")
samples_table <- read.table(samples_file,header=T)
samples_list <- unique(samples_table$sample)

log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# Now sum every three rows
print(dim(counts_with_reps))
master <- as.data.frame(sapply(seq(2,dim(counts_with_reps)[2],by=3),
                               function(i) rowSums(counts_with_reps[,i:(i+2)])))

# We removed the first column, "genes", from the last function, and we want it as our rownames

rownames(master) <- counts_with_reps$gene
print("this many samples")
print(samples_list)
print(length(samples_list))
print("Master is:\t")
print(dim(master))
print(colnames(master))

master$gene <- rownames(master)
colnames(master) <- c(samples_list,"gene")
library(tidyverse)
master <- master %>% select(c("gene",samples_list))

write.table(master,out_filename,quote=F,sep="\t",row.names = F)

# Label "master" with list of sample names \

# names(master) <- list.of.sample.names

# Split the experiments 
#transgenic.exp <- master[,1:47]
#timecourse.exp <- master[,48:135]

# Write output files for each exp, and you're done
