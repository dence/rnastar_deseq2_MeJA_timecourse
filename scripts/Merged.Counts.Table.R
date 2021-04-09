# Make a Merged Counts Table

# List of inputs you need:

# spike_in.counts <- # need input table with ERCC info by lane
# counts_with_reps <- # need input table with gene info by lane
# list.of.sample.names <- # need input list of sample names as character vector with length 135

# First, rbind the tables together

all.counts <- rbind(spike_in.counts,counts_with_reps)

# Now sum every three rows

master <- as.data.frame(sapply(seq(2,406,by=3),function(i) rowSums(all.counts[,i:(i+2)])))

# We removed the first column, "genes", from the last function, and we want it as our rownames

rownames(master) <- all.counts$gene

# Label "master" with list of sample names \

# names(master) <- list.of.sample.names

# Split the experiments 

transgenic.exp <- master[,1:47]
timecourse.exp <- master[,48:135]

# Write output files for each exp, and you're done
