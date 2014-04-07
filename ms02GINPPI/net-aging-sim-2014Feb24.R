# 2014 Feb 17, change name "20131221.DIPandGIN.sim.aging_v2.R" to "net-aging-sim-2014Feb17.R"

# 2013 Dec 20, merge DIP PPI and Genetic Inxt Net -> Multi-net approach
rm(list=ls())

require('flexsurv')
#source("/Users/hongqin/lib/R/lifespan.r")
source("lifespan.r")
source("network.r")

#setwd("~/projects/0.network.aging.prj/0.ppi.reliability.simulation")
list.files(path='data', )
setwd("~/projects/0.ginppi.reliability.simulation/ms02GINPPI")

debug = 0; 

#essential gene info
essenTb = read.csv('SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))

# remove self-intxns

#for ( i in 1:100) #
for( i in 3:100 ){ ##start at 3 after power outage, 20140407
  path = paste('dipgin.ms02.output/', i, sep='')
  ms02file = paste('ms02_', i, ".tab", sep='')
  infile = paste( path, '/', ms02file, sep=""); print(infile)
  pairs = read.table( infile, header=T, sep="\t", colClass = c("character", "character", NA) )
  print(head(pairs))
  if(debug==9) {     pairs = pairs[1:1000,]  }
  pairs = pairs[ pairs$id1 != pairs$id2, ]  

  # How do the two data set overlap? DIP seems to contain some questionable orfs
  uniq.orf.from.pairs = unique(c(pairs$id1, pairs$id2)) #5496 ORF
  matches = uniq.orf.from.pairs %in% unique(essenTb$orf)
  table(matches)
  #unmatchedORF = uniq.orf.from.pairs[! matches]

  #matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='essential'])
  #table(matches)  
  #matches = uniq.orf.from.pairs %in% unique(essenTb$orf[essenTb$essenflag=='nonessential'])
  #table(matches)
  #remove unmatched orfs from pairs
  #pairs$Removeflag = ifelse( pairs$id1 %in%unmatchedORF | pairs$id2 %in%unmatchedORF, T,F   )
  #table(pairs$Removeflag)
  #pairs = pairs[! pairs$Removeflag, ]
  #table(pairs$Removeflag)
  #pairs = pairs[,1:2]  ##This set of pairs is read for analysis
  
  # label essential nodes, remove nonesse-nonessen pairs
  essentialORFs = essenTb$orf[essenTb$essenflag=='essential']
  pairs$essen1 = pairs$id1 %in% essentialORFs
  pairs$essen2 = pairs$id2 %in% essentialORFs
  #remove nonessen <-> nonessen intxn because they do not affect aging. 
  pairs$remove = ifelse( pairs$essen1==F & pairs$essen2==F, T, F  )
  pairs= pairs[! pairs$remove, ]  
  # 34052 for one ms02 network
  
  #how many essen <--> essen intxn? 
  pairs$inxnEE = pairs$essen1 & pairs$essen2
  table(pairs$inxnEE)
  #How many essen genes? 
  tmp = essentialORFs %in% unique(c(pairs$id1, pairs$id2)) 
  table(tmp)
  essentialORFsPPI = essentialORFs[tmp] #?????
  
  #get connectivities per node
  degreeTb = data.frame( table(c(pairs[,1], pairs[,2])))
  summary(degreeTb)
  #median degree =5, mean=12
  #for one ms02, media =6, mean=13.68, so orginal network is power-law like, skew at two ends. 
  degreeTb$ORF = as.character( degreeTb[,1])
  #hist(log2(degreeTb$Freq), breaks=30)
  
  degreeCutoff = 4; #######!!!!!! degree=5 is the median
  tmp = essentialORFsPPI %in% degreeTb$ORF[degreeTb$Freq>degreeCutoff]
  GooddEssentialORFsPPI = essentialORFsPPI[tmp]
    
  if(debug >= 5){GooddEssentialORFsPPI = GooddEssentialORFsPPI[1:100]  }
  
  lambda_v = 1/c(100, 50, 25 ) #2014Feb 25, fix double inverse bug in lambda
  p_v = seq(0.6, 1.0, by=0.1)  ; #the chance that each gene interaction is active at t=0
  sim_names = c( "degreeCutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
  sim       = t( c(NA,     NA,   NA,       NA,       NA,        NA,  NA,   NA,      NA))
  sim = data.frame(sim)
  names(sim) = sim_names

  full_age_dir = paste(path, '/', 'popages', sep='')
  system(paste('mkdir ', full_age_dir ))
  
  for(lambda in lambda_v) {  
    for( p in p_v) {  # p=0.9, #for debug
      #popSize = 500 #too small pop size and too small p can lead to very few living individuals
      popSize = 2000 #20140317 Monday
      popAges = numeric(popSize)
      time1 = date()
      j=1; count = 0; 
      while ((j <= popSize) && ( count < popSize*30)) {
        count = count + 1;      
        print(paste("count=",count))
        currentNetworkAge = single_network_failure(lambda, p, pairs, GooddEssentialORFsPPI)
        #single_network_failure = function(lambda, p, pairs, runningORFs) {   
        if (currentNetworkAge > 0) {
          popAges[j] = currentNetworkAge      
          j = j+1
        } 
      }# end of j while-loop, population loop
            
      
      timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
      age.file.name=paste("cutoff", degreeCutoff, "p", p, "lambda", lambda, 'popsize',popSize, "time",timestamp, "txt", sep="." )
      full_age_file = paste( full_age_dir,'/', age.file.name, sep='')
      
      write.csv( popAges, full_age_file, row.names=F)
      
      # s.tb = calculate.s ( popAges )      #plot( s.tb$s ~ s.tb$t, type='l', log='x' )       
      # lifespanGomp = flexsurvreg(formula = Surv(popAges) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
      # lifespanWeib = flexsurvreg(formula = Surv(popAges) ~ 1, dist = 'weibull')  
      # c(lifespanWeib$AIC, lifespanGomp$AIC, lifespanWeib$AIC - lifespanGomp$AIC )
      # sOject = Surv(popAges)
      # sim_names = c( "cutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
      # sim = rbind(sim, c( degreeCutoff, p, lambda, mean(popAges), median(popAges), 
      #                    lifespanGomp$res[2,1], lifespanGomp$res[1,1], lifespanGomp$AIC, lifespanWeib$AIC))
      
    } # end of p-loop  
    #timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
    #write.csv(sim, file= paste(currentwkdir, "sce_sim_", timestamp, ".csv", sep=""), row.names=F)
  } #end of lambda loop
  
}  
