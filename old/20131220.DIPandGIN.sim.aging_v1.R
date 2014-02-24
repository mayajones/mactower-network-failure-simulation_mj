# todo, 
# permutation effect on aging? 
# lambda ~ 1/connectivity of nodes

# 2013 Dec 20, merge DIP PPI and Genetic Inxt Net -> Multi-net approach
rm(list=ls())

#require(date)
require('flexsurv')
#source("/Users/hongqin/lib/R/lifespan.r")
source("lifespan.r")

#setwd("~/projects/0.network.aging.prj/0.ppi.reliability.simulation")
list.files(path='data', )

debug = 0; 

#yeast PPI
#pairs = read.csv('data/pairs.csv', colClasses=c('character','character'))
#this yeast ppi dataset is consistent with Taiwan group's report.
dip = read.csv("data/yeastDIP.csv")
pairsPPI = dip[,c(1,2)]
pairsPPI$ORF1 = as.character(pairsPPI$ORF1)
pairsPPI$ORF2 = as.character(pairsPPI$ORF2)
pairsPPI = pairsPPI[ pairsPPI$ORF1 != pairsPPI$ORF2, ]#remove self-intxns

# yeast genetic network
#gPairs = read.csv("data/sgadata_costanzo2009_stringentCutoff_101120.csv", header=F)
gPairs = read.csv("data/sgadata_costanzo2009_lenientCutoff_101120.csv", header=F)
names(gPairs) = c("ORF1", "Name1", "ORF2", "Name2", NA, NA, NA)
gPairs$ORF1 = as.character( gPairs$ORF1 )
gPairs$ORF2 = as.character( gPairs$ORF2 )

#merge PPI and GIN
pairs = rbind(pairsPPI[,c("ORF1","ORF2")], gPairs[,c("ORF1","ORF2")])

# remove self-intxns
pairs = pairs[ pairs$ORF1 != pairs$ORF2, ]
# 96851 pairs for DIP+GIN.strigentCutoff
# 786118 for DIP+GIN.lenientCutoff

#essential gene info
essenTb = read.csv('data/SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))

#######################
# How do the two data set overlap? DIP seems to contain some questionable orfs
uniq.orf.from.pairs = unique(c(pairs$ORF1, pairs$ORF2)) #4207 ORF
matches = uniq.orf.from.pairs %in% unique(essenTb$orf)
table(matches)
#FALSE  TRUE 
# 720  5507

unmatchedORF = uniq.orf.from.pairs[! matches]

matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='essential'])
table(matches)
#FALSE  TRUE 
#5171   1056  
#This is a good coverage

matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='nonessential'])
table(matches)
# FALSE  TRUE 
# 1839  4388  #this is amazingly consistent with Taiwan group's report.

#remove unmatched orfs from pairs
pairs$Removeflag = ifelse( pairs$ORF1 %in%unmatchedORF | pairs$ORF2 %in%unmatchedORF, T,F   )
table(pairs$Removeflag)
# FALSE  TRUE 
# 572221  213897

#So, the updated DIP has         19770 intxn
#So, the DIP+GIN.lenior lead to 572221 intxns

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
pairs= pairs[! pairs$remove,1:4 ]  
#Now, only 31394 intxn left

#how many essen <--> essen intxn? 
pairs$inxnEE = pairs$essen1 & pairs$essen2
table(pairs$inxnEE)
# FALSE  TRUE 
# 28222  3172
#So, EE intxn should only occur in PPI

#How many essen genes? 
tmp = essentialORFs %in% unique(c(pairs$ORF1, pairs$ORF2)) 
table(tmp)
#FALSE  TRUE 
#110   1048  

essentialORFsPPI = essentialORFs[tmp] #?????

#get connectivities per node
degreeTb = data.frame( table(c(pairs[,1], pairs[,2])))
summary(degreeTb)
#median degree =5, mean=12

degreeTb$ORF = as.character( degreeTb[,1])
hist(log2(degreeTb$Freq), breaks=30)

degreeCutoff = 4; #######!!!!!! degree=5 is the median
tmp = essentialORFsPPI %in% degreeTb$ORF[degreeTb$Freq>degreeCutoff]
GooddEssentialORFsPPI = essentialORFsPPI[tmp]

###########################

if(debug >= 5){GooddEssentialORFsPPI = GooddEssentialORFsPPI[1:100]  }

###########################
# simulate aging
# -> exponential age to all pairs
# -> maximal age for each essential gene
# -> minimal age for all essential modules

set.seed(2013)

lambda_v = 1 / 3^seq(3,5)
#lambda_v = 1 / s

# Notice that degree*p > 1 if the cells are to be alive. >2 should sufficient redundancy. 
p_v = seq(0.4, 1.0, by=0.1)  ; #the chance that each gene interaction is active at t=0
#p_v = c(0.7, 0.7, 0.7,  0.8, 0.8, 0.8, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0)  ; #the chance that each gene interaction is active at t=0

sim_names = c( "degreeCutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
sim       = t( c(NA,     NA,   NA,       NA,       NA,        NA,  NA,   NA,      NA))
sim = data.frame(sim)
names(sim) = sim_names


################################## network simulations
# for debug, lambda = 1/10
# lambda = 1/35; p=0.9

for(lambda in lambda_v) {
  # p=0.9, #for debug
  for( p in p_v) {  
    popSize = 100 #too small pop size and too small p can lead to very few living individuals
    popAges = numeric(popSize)
    time1 = date()
    j=1; count = 0; 
    while ((j <= popSize) && ( count < popSize*100)) {
      count = count + 1; 
      runningORFs = GooddEssentialORFsPPI  
      
      #I could add stochasticity into pairs here.  
      pairs$active = runif(length(pairs[,1]))
      # tmp = pairs$active > 1-p
      # table(tmp) / length(tmp)  ; #double-check, very good. 
      
      ModuleTb = data.frame(runningORFs)
      pairs$age = rexp( length(pairs[,1]), rate=lambda )  #exponential ages for pairs
      pairs$age = ifelse(pairs$active > (1-p), pairs$age, NA )
      
      #loop every essential genes to identify the module age
      for (i in 1:length(runningORFs)) {
        myORF = runningORFs[i]
        pos1 = grep(myORF, pairs$ORF1)
        pos2 = grep(myORF, pairs$ORF2)
        if( length( c(pos1,pos2))>=1 ) {
          ModuleTb$age.m[i] = max( pairs$age[c(pos1,pos2)], na.rm=T )   #maximal intxn age -> module age
        } else {
          ModuleTb$age.m[i] = NA; 
        }
      }
      #head(ModuleTb); 
      summary(ModuleTb)
      ModuleTb$age.m[ ModuleTb$age.m== -Inf] = 0; #dead births occur when links are not active
      currentNetworkAge = min(ModuleTb$age.m)
      if (currentNetworkAge > 0) {
        popAges[j] = currentNetworkAge      
        j = j+1
      } 
    }# end of j while-loop, population loop
    
    time2 = date()
    hist(popAges)
    summary(popAges)
    popAges = popAges[popAges>0]; #remove dead-births, which can occur when p is low
    
    #time1; time2; 
    #s.tb = calculate.s ( popAges )
    #plot( s.tb$s ~ s.tb$t ) 
    #plot( s.tb$s ~ s.tb$t, type='l', log='x' ) 
        
    lifespanGomp = flexsurvreg(formula = Surv(popAges) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
    lifespanWeib = flexsurvreg(formula = Surv(popAges) ~ 1, dist = 'weibull')  
    c(lifespanWeib$AIC, lifespanGomp$AIC, lifespanWeib$AIC - lifespanGomp$AIC )
    sOject = Surv(popAges)
    
    age.file.name=paste("cutoff", degreeCutoff, "p", p, "lambda", lambda, as.character(date()), "txt", sep="." )
    write.csv( popAges, paste("simulated.ages/",age.file.name, sep=""), row.names=F)
    
    #sim_names = c( "cutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
    sim = rbind(sim, c( degreeCutoff, p, lambda, mean(popAges), median(popAges), 
                        lifespanGomp$res[2,1], lifespanGomp$res[1,1], lifespanGomp$AIC, lifespanWeib$AIC))
    
  } # end of p-loop
  
  timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
  write.csv(sim, file= paste("sceGINPPIaging_sim_", timestamp, ".csv", sep=""), row.names=F)
  
} #end of lambda loop

write.csv(sim, file="sceGINPPIaging_sim_2013121_end.csv", row.names=F)

#quite('yes')
