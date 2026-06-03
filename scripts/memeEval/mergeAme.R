#This script forms the overlap of enriched TFBS from real- and shuffled controls

###
###
###Libraries
set.seed(1234)
library(data.table)

###
###
###Main

##
##Read input and merge
realCon <- fread("../../data/memeEval/ame/caABC_pancreas07_mergedRepli_digestiveHpa25/ameOutHocomoco/ame.tsv",header = TRUE, data.table = FALSE, select = c(3,8,15,17))
shufCon <- fread("../../data/memeEval/ame/caABC_pancreas07_mergedRepli_digestiveHpa25/ameOutHocomocoShuf/ame.tsv",header = TRUE, data.table = FALSE, select = c(3,8,15,17))
colnames(realCon) <- c("motif_ID", "E_value_realCon", "%TP_realCon", "%FP_realCon")
colnames(shufCon) <- c("motif_ID", "E_value_shufCon", "%TP_shufCon", "%FP_shufCon")
mergedData <- merge(realCon,shufCon)
mergedData <- mergedData[order(mergedData$E_value_realCon),]

##
##Write output
write.table(x = mergedData, file = "../../data/memeEval/ame/caABC_pancreas07_mergedRepli_digestiveHpa25/mergedData.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
