# todo, 
# permutation effect on aging? 
# lambda ~ 1/connectivity of nodes

rm(list=ls())

require('flexsurv')
#source("/Users/hongqin/lib/R/lifespan.r")
source("lifespan.r")

#setwd("~/projects/0.network.aging.prj/0.ppi.reliability.simulation")
list.files(path='data', )

debug = 0; 

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
#3506   972  #this is amazingly consistent with the Taiwan group's report.

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


degreeCutoff = 5; #######!!!!!!
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

lambda_v = 1 / 3^seq(2,5)
#lambda_v = 1 / s

#p_v = seq(0.1, 1.0, by=0.1)  ; #the chance that each gene interaction is active at t=0
p_v = c(0.7, 0.7, 0.7,  0.8, 0.8, 0.8, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0)  ; #the chance that each gene interaction is active at t=0

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
    for( j in 1:popSize) {
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
      popAges[j] = currentNetworkAge      
    }# end of j loop, population loop
    
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
    
    #sim_names = c( "cutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
    sim = rbind(sim, c( degreeCutoff, p, lambda, mean(popAges), median(popAges), 
                        lifespanGomp$res[1,1], lifespanGomp$res[2,1], lifespanGomp$AIC, lifespanWeib$AIC))
    
  } # end of p-loop
  
  write.csv(sim, file="scePPIaging_sim_20131217_v2_cutoff5_tmp.csv", row.names=F)
  
} #end of lambda loop

#conditions$avgLS[r] = mean(lifespansTemp)
#conditions$medianLS[r] = median(lifespansTemp)
#conditions$gompShape[r] = lifespanGomp$res[1,1]
#conditions$gompRate[r] = lifespanGomp$res[2,1]
#conditions$gompLogLik[r] = lifespanGomp$loglik
#conditions$gompAIC[r] = lifespanGomp$AIC
#conditions$weibShape[r] = lifespanWeib$res[1,1]
#conditions$weibScale[r] = lifespanWeib$res[2,1]
#conditions$weibLogLik[r] = lifespanWeib$loglik
#conditions$weibAIC[r] = lifespanWeib$AIC    

write.csv(sim, file="scePPIaging_sim_20131217_v2_cutoff5_end.csv", row.names=F)

#quite('yes')
