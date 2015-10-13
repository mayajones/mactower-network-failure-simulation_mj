# 20151012 Plan 
# input a gene network and essential gene list. 
# remove self-interactions
# assign smaller number to essential nodes, and large numbers to non-essential nodes
# identify essential and non-essential nodes in gene networks
# remove non-essential gene -noessential gene interactions.

rm(list=ls())
source("lifespan.r")
source("network.r")

setwd("~/github/mactower-network-failure-simulation")
list.files(path='data', )
debug = 0; 

#essential gene info
 #essenTb = read.csv('SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))
 #20151012, fixed a bug due to comma in orf. 
essenTb = read.csv("SummaryRegressionHetHom2015Oct12.csv", colClasses=rep('character', 9))
length(unique(essenTb$orf)) #5772,

  infile = "data/merged_PPIGIN_2014Jan20.tab";
  pairs = read.table( infile,  header=F, sep="\t", colClass = c("character", "character") )
  names(pairs) = c("id1",'id2')
  print(head(pairs))
  if(debug==9) {     pairs = pairs[1:1000,]  }
  pairs = pairs[ pairs$id1 != pairs$id2, ]  

  # How do the two data set overlap? DIP seems to contain some questionable orfs
  uniq.orf.from.pairs = unique(c(pairs$id1, pairs$id2)) #5496 ORF
  matches = uniq.orf.from.pairs %in% unique(essenTb$orf)
  table(matches)
  #unmatchedORF = uniq.orf.from.pairs[! matches]

  # Generate geneNOs (numeric IDs) using the essential, nonessential gene file. 
  essentialORFs = essenTb$orf[essenTb$essenflag == 'essential']
  essentialORFs = essentialORFs[order(essentialORFs)]
  nonessentialORFs = essenTb$orf[essenTb$essenflag != 'essential']
  mergedORFs = c(essentialORFs, nonessentialORFs)
  mergedNOs = 1: length(mergedORFs)
  names(mergedNOs) = mergedORFs   # This is the numbers assigned to ORFs/Nodes !!!
  mergedNOs[c('YAL016W', "YPR116W")] #check, passed 20151012Mon
  essenTb$geneNO = mergedNOs[essenTb$orf]
  essenTb$essenflagNumeric = ifelse( essenTb$essenflag=='essential', 1, 0)
  rownames(essenTb) = essenTb$orf
  head(essenTb)
  
  essenTb2 = essenTb[, c("geneNO", "essenflagNumeric")]
  essenLookupTb = essenTb2$essenflagNumeric[order(essenTb2$geneNO)]
  #manual check
  essenLookupTb[c(1152,1153)]
  
  ### convert pairs ORF into geneNOs
  essenTb[pairs$id1[1:10], 'geneNO']
  pairs$No1 = essenTb$geneNO[match(pairs$id1, essenTb$orf)]
  pairs$No2 = essenTb$geneNO[match(pairs$id2, essenTb$orf)]
  #manually check a few
  pairs[100, ]
  # id1     id2 No1  No2
  # 100 YDL043C YGL049C 131 2623
  essenTb["YDL043C", ]
  essenTb["YGL049C", ]
  
  ##output
  write.csv(essenTb,"SummaryRegressionHetHomFactorized2015Oct13.csv")
  write.csv(pairs, "data/merged_PPIGIN_Factorized2015Oct13.csv")
  
