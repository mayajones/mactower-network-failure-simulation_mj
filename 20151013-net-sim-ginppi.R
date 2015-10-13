
##########################################
#2015Oct13, use numeric lookup table for essential genes.
#2014March3, redo ginppi simulation witht same parameter for ms02networks. 
# 2014 Feb 17, change name "20131221.DIPandGIN.sim.aging_v2.R" to "net-aging-sim-2014Feb17.R"
# 2013 Dec 20, merge DIP PPI and Genetic Inxt Net -> Multi-net approach
#rm(list=ls())

single_network_failure_v2 = function(lambda1, lambda2=lambda1/10, threshold=4, p, pairs, essenLookupTb ) {
  # single network failure simulation, 20151013Tue
  # lambda1: First exponential constant failure rate for edges with degree > threshold
  # lambda2: Second exponential constant failure rate for edges with degree <= threshold
  # threshold: degree threshold for lambda1 and lambda2
  # pairs: network in pairwide format, using numeric NOs 20151013
  # essenLookupTb: lookup table for essential and nonessential genes, numeric values 
  ## for debug:   lambda1 = 1/50; lambda2= lambda1/10; threshold=4; p=0.8
  
  inpairs = pairs[,3:4] #bookkeeping  
  names(inpairs) = c('No1','No2')
  
  #get connectivities per node
  degreeTb = data.frame( table(c(inpairs$No1, inpairs$No2)))
  names(degreeTb) = c("No", "degree")
  degreeTb$moduleAge = NA;
  
  for( i in 1:length(degreeTb[,1])){
    if ( essenLookupTb[ degreeTb$No[i] ]) { #essential node
      lambda = ifelse( degreeTb$degree[i] >= threshold, lambda1, lambda2)
      age = rexp( degreeTb$degree[i], rate=lambda ) #exponential age
      if(degreeTb$degree[i] >= threshold){
        active = runif(degreeTb$degree[i])  #uniform interaction stochasticity
        active = ifelse( active<=p, 1, NA  ) #pick active interactions
        if( sum(active, na.rm=T) > 0 ){ #there should be at least 1 active intxn
          age = age * active # only active interactions for modular age estimation
          degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
        } else {# when no active intxn is available 
          degreeTb$moduleAge[i] = 0; #this module is born dead.
        }
      } else { # for degree < threshold, no stochasticity is applied. 
        degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
      }
    } else {# non-essential node
      degreeTb$moduleAge[i] = NA 
    }
  }
  
  summary(degreeTb)
  currentNetworkAge = min(degreeTb$moduleAge, na.rm=T)
}

#source("network.r")

# R -f file --args lambda1 lambda2 degreeThreshold p popSize
# R -f 20151013-net-sim-ginppi.R --args 0.002 0.0002 4 0.9 5
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
lambda1 = as.numeric(args[1]); lambda1;
lambda2 = as.numeric(args[2]); lambda2;
degreeThreshold = as.integer((args[3])); degreeThreshold
p = as.numeric(args[4]); p;
popSize = as.numeric(args[5]); popSize;

list.files(path='data', )
debug = 0; 

#essential gene info
essenTb = read.csv("data/SummaryRegressionHetHomFactorized2015Oct13.csv", colClasses=rep('character', 9))
essenLookupTb = read.csv("data/essntialGeneLookupTable_20151013.csv", row.names=1)
essenLookupTb = essenLookupTb[,1]

infile = "data/merged_PPIGIN_Factorized2015Oct13.csv"
pairs = read.csv(infile)
  names(pairs) = c("id1",'id2', "No1", "No2")
  print(head(pairs))
  if(debug==9) {     pairs = pairs[1:1000,]  }
  pairs = pairs[ pairs$No1 != pairs$No2, ]  

  # label essential nodes, remove nonesse-nonessen pairs
  pairs$essen1 = essenLookupTb[pairs$No1]
  pairs$essen2 = essenLookupTb[pairs$No2]
  #remove nonessen <-> nonessen intxn because they do not affect aging. 
  pairs$remove = ifelse( pairs$essen1==F & pairs$essen2==F, T, F  )
  pairs= pairs[! pairs$remove, ]  
  # 31394 for one ms02 network
  
  #get connectivities per node
  degreeTb = data.frame( table(c(pairs$No1, pairs$No2)))
  summary(degreeTb); 
  degreeTb[1:10,]
  #median degree =5, mean=12
  #for one ms02, media =6, mean=13.68, so orginal network is power-law like, skew at two ends. 

  full_age_dir = "ori.ginppit.2015Oct"  
  
      popAges = numeric(popSize)
      time1 = date()
      j=1; count = 0; 
      while ((j <= popSize) && ( count < popSize*30)) {
        count = count + 1;      
        print(paste("count=",count))
        currentNetworkAge = single_network_failure_v2(lambda1, lambda2, degreeCutoff, p, pairs, essenLookupTb)
        if (currentNetworkAge > 0) {
          popAges[j] = currentNetworkAge      
          j = j+1
        } 
      }# end of j while-loop, population loop
            
      
      timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
      age.file.name=paste("threshold", degreeCutoff, "p", p, "lambda1", lambda1, 
                          "lambda2", lambda2,'popsize',popSize, "time",timestamp, "txt", sep="." )
      full_age_file = paste( full_age_dir,'/', age.file.name, sep='')
      
      write.csv( popAges, full_age_file, row.names=F)
      
