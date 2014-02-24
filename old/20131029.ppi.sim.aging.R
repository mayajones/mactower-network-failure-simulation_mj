rm(list=ls())

list.files(path='data', )

#yeast PPI
pairs = read.csv('data/pairs.csv', colClasses=c('character','character'))
#this yeast ppi dataset is consistent with Taiwan group's report.

#essential gene info
essenTb = read.csv('data/SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))

#######################
# How do the two data set overlap? DIP seems to contain some questionable orfs
uniq.orf.from.pairs = unique(c(pairs$ORF1, pairs$ORF2)) 
matches = uniq.orf.from.pairs %in% unique(essenTb$orf)
table(matches)
#FALSE  TRUE 
# 261  4217
# YDL026W dubious, YAR009C transposable element
dubiousORF = uniq.orf.from.pairs[! matches]

matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='essential'])
table(matches)
#FALSE  TRUE 
#3506   972  #this is amazingly consistent with Taiwan group's report.

matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='nonessential'])
table(matches)
# FALSE  TRUE 
# 1287  3191  #this is amazingly consistent with Taiwan group's report.

#remove dubious orfs from PPI
pairs$Removeflag = ifelse( pairs$ORF1 %in%dubiousORF | pairs$ORF2 %in%dubiousORF, T,F   )
table(pairs$Removeflag)
#FALSE  TRUE 
#12180  1487 

pairs = pairs[! pairs$Removeflag, ]
table(pairs$Removeflag)
pairs = pairs[,1:2]  ##This set of pairs is read for analysis

###############################
# label essential nodes, remove nonesse-nonessen pairs
essentialORFs = essenTb$orf[essenTb$essenflag=='essential']
pairs$essen1 = pairs$ORF1 %in% essentialORFs
pairs$essen2 = pairs$ORF2 %in% essentialORFs

head(pairs)

#remove nonessen <-> nonessen intxn because they do not affect aging. 
pairs$remove = ifelse( pairs$essen1==F & pairs$essen2==F, T, F  )
pairs= pairs[! pairs$remove,1:4 ]  #only 6279 intxn left

#remove self-intxn, just to make sure
pairs = pairs[ pairs$ORF1 != pairs$ORF2, ]

#how many essen <--> essen intxn? 
pairs$inxnEE = pairs$essen1 & pairs$essen2
table(pairs$inxnEE)
# FALSE  TRUE 
# 4521  1758 

#How many essen genes? 
tmp = essentialORFs %in% unique(c(pairs$ORF1, pairs$ORF2)) 
table(tmp)
#FALSE  TRUE 
#194   958  #So, intxn among essential genes appear to be suppressed. 

essentialORFsPPI = essentialORFs[tmp]

#get connectivities per node
degreeTb = data.frame( table(c(pairs[,1], pairs[,2])))
summary(degreeTb)
degreeTb$ORF = as.character( degreeTb[,1])

degreeCutoff = 4; #######!!!!!!
tmp = essentialORFsPPI %in% degreeTb$ORF[degreeTb$Freq>degreeCutoff]
GooddEssentialORFsPPI = essentialORFsPPI[tmp]


###########################
# simulate aging
# -> exponential age to all pairs
# -> maximal age for each essential gene
# -> minimal age for all essential genes

set.seed(2013)

lambda = 1/3500

# runningORFs = essentialORFsPPI
runningORFs = GooddEssentialORFsPPI  

popSize = 100
popAges = numeric(popSize)
time1 = date()
for( j in 1:popSize) {
  ModuleTb = data.frame(runningORFs)
  pairs$age = rexp( length(pairs[,1]), rate=lambda )  #exponential ages for pairs
  
  #I could add stochasticity into pairs here.  
  
  for (i in 1:length(runningORFs)) {
  myORF = runningORFs[i]
  pos1 = grep(myORF, pairs$ORF1)
  pos2 = grep(myORF, pairs$ORF2)
  ModuleTb$age.m[i] = max( pairs$age[c(pos1,pos2)] )  
 }
 head(ModuleTb); 
  summary(ModuleTb)
  currentNetworkAge = min(ModuleTb$age.m)
  popAges[j] = currentNetworkAge
}
time2 = date()
hist(popAges)
summary(popAges)

time1; time2; 
source("/Users/hongqin/lib/R/lifespan.r")
s.tb = calculate.s ( popAges )
plot( s.tb$s ~ s.tb$t ) 
# This is exponential, perhaps not surprisingly, because power-law elevate the role of the weakest link. 
# So, if the weakest links die more slowly, maybe we can see sigmoidal shapes. 
# Maybe I should remove the single-linked essential genes. In essence, assumming they die very slow. 
# maybe I should also use house-keeping genes

plot( s.tb$s ~ s.tb$t, type='l', log='x' ) 

# todo, 
# permutation effect on aging? 
# lambda ~ 1/connectivity of nodes

# I should try to normalized the survival curves for a different prespective



#quite('yes')
