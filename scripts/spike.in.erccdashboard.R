# Spike in Normalization Script

# erccdashboard

library(erccdashboard)

# Import spike_in.counts.tsv file
all <- read.delim("/blue/peter/mallory.morgan/rnastar_deseq2_MeJA_timecourse/results/counts/spike_in.counts.tsv")

# Remove "sample_" and replace with ""
names(all) = gsub(pattern = "sample_*", replacement = "", x = names(all))

# Subset files into two experiments ###CHANGE THIS IF INTERMIXED###
trg = subset(all, select = c(1:142)) # transgenics (47 samples x 3 lanes per sample)
tc = subset(all, select = c(1,143:407)) # timecourse (88 samples x 3 lanes per sample)

# Rename with trt and controls for each experiment

# transgenics
vc = c("277-L1","277-L2","277-L3","278-L1","278-L2","278-L3","279-L1","279-L2","279-L3","280-L1","280-L2","280-L3") # pWVK305 samples (vector controls)
names(trg) <- ifelse(names(trg) %in% vc ,paste("CTL", colnames(trg), sep = "_"),paste("TRT", colnames(trg), sep = "_"))
names(trg)[1] <- "genes"

# timecourse
controls = c("M2-L1","M2-L2","M2-L3","M21-L1","M21-L2","M21-L3","M23-L1","M23-L2","M23-L3","M34-L1","M34-L2","M34-L3",
             "M36-L1","M36-L2","M36-L3","M49-L1","M49-L2","M49-LL3","M54-L1","M54-L2","M54-L3","M59-L1","M59-L2",
             "M59-L3","M71-L1","M71-L2","M71-L3","M76-L1","M76-L2","M76-L3","M8-L1","M8-L2","M8-L3","M80-L1","M80-L2","M80-L3",
             "M86-L1","M86-L2","M86-L3","M90-L1","M90-L2","M90-L3","M94-L1","M94-L2","M94-L3") # All control days are compiled
names(tc) <- ifelse(names(tc) %in% controls, paste("CTL", colnames(tc), sep = "_"), paste("TRT", colnames(test), sep = "_"))
names(tc)[1] <- "genes"

# transgenic analysis

exDat <- runDashboard(datType = "count", isNorm = FALSE, exTable = trg,
                      filenameRoot = "transgenic", sample1Name = "TRT", sample2Name = "CTL",
                      erccmix = "Single", erccdilution = 1/1000, spikeVol = 1, totalRNAmass = 0.05,
                      choseFDR = 0.05)
# timecourse
exDat <- runDashboard(datType = "count", isNorm = FALSE, exTable = tc,
                      filenameRoot = "timecourse", sample1Name = "TRT", sample2Name = "CTL",
                      erccmix = "Single", erccdilution = 1/1000, spikeVol = 1, totalRNAmass = 0.05,
                      choseFDR = 0.05)

# datType = "count" # "count" for RNA-Seq data, "array" for microarray data
# isNorm = FALSE # flag to indicate if input expression measures are already
# normalized, default is FALSE
# exTable = MET.CTL.countDat # the expression measure table
# filenameRoot = "RatTox" # user defined filename prefix for results files
# sample1Name = "MET" # name for sample 1 in the experiment
# sample2Name = "CTL" # name for sample 2 in the experiment
# erccmix = "RatioPair" # name of ERCC mixture design, "RatioPair" is default
# erccdilution = 1/100 # dilution factor used for Ambion spike-in mixtures
# spikeVol = 1 # volume (in microliters) of diluted spike-in mixture added to
# total RNA mass
# totalRNAmass = 0.500 # mass (in micrograms) of total RNA
# choseFDR = 0.05 # user defined false discovery rate (FDR), default is 0.05
