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
vc = c("277", "278", "279", "280") # pWVK305 samples (vector controls)
names(trg) <- ifelse(names(trg) %in% vc ,paste("CTL", colnames(trg), sep = "_"),paste("TRT", colnames(trg), sep = "_"))
names(trg)[1] <- "genes"

# timecourse
controls = c("M2","M21","M23","M34","M36","M49","M54",
             "M59","M71","M76","M8","M80","M86","M90","M94") # All control days are compiled
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
