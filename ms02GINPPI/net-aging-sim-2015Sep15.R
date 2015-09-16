#20150915Mon add argument for lambda, p, popSize

# 2014 Feb 17, change name "20131221.DIPandGIN.sim.aging_v2.R" to "net-aging-sim-2014Feb17.R"

# 2013 Dec 20, merge DIP PPI and Genetic Inxt Net -> Multi-net approach
rm(list=ls())

source("lifespan.r")
source("network.r")

#for debug
#  start =1; end = 1; lambda= 1/100; p = 0.95; popSize=5;


# R -f file --args start end lambda p popSize
#  R -f net-aging-sim-2014Sep15.R --args 1 1 0.002 0.95 10
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))
start = as.integer(args[1]); start; 
end = as.integer(args[2]); end; 
lambda = as.numeric(args[3]); lambda;
p = as.numeric(args[4]); p;
popSize = as.numeric(args[5]); popSize;


myhost = 'greenfield'  # 'byte' 'blacklight' 'mactower'
#myhost = 'byte'  # 'byte' 'blacklight' 'mactower'
#myhost = 'helen';

mydir = "/crucible/mc48o9p/hqin2/mactower-network-failure-simulation-master/ms02GINPPI"
if (myhost == 'byte') {  mydir = "/Users/hqin/github/mactower-network-failure-simulation/ms02GINPPI"
} else if (myhost == 'helen') { mydir = "/Users/hqin/github/mactower-network-failure-simulation/ms02GINPPI";  
} 

mydir
setwd(mydir)

debug = 0; 

#essential gene info
essenTb = read.csv('SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))


for ( i in start:end) {
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
  degreeTb$ORF = as.character( degreeTb[,1])
  
  degreeCutoff = 4; #######!!!!!! degree=5 is the median
  tmp = essentialORFsPPI %in% degreeTb$ORF[degreeTb$Freq>degreeCutoff]
  GooddEssentialORFsPPI = essentialORFsPPI[tmp]
    
  if(debug >= 5){GooddEssentialORFsPPI = GooddEssentialORFsPPI[1:100]  }
  
  sim_names = c( "degreeCutoff","p", "lambda", "meanLS", "medianLS", "R","G", "GompAIC", "WeibAIC")
  sim       = t( c(NA,     NA,   NA,       NA,       NA,        NA,  NA,   NA,      NA))
  sim = data.frame(sim)
  names(sim) = sim_names

  full_age_dir = paste(path, '/', 'popages', sep='')
  system(paste('mkdir ', full_age_dir ))
  
      popAges = numeric(popSize)
      time1 = date()
      j=1; count = 0; 
      while ((j <= popSize) && ( count < popSize*30)) {
        count = count + 1;      
        print(paste("count=",count))
        currentNetworkAge = single_network_failure(lambda, p, pairs, GooddEssentialORFsPPI)
        if (currentNetworkAge > 0) {
          popAges[j] = currentNetworkAge      
          j = j+1
        } 
      }# end of j while-loop, population loop
            
      
      timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
      age.file.name=paste("cutoff", degreeCutoff, "p", p, "lambda", lambda, 'popsize',popSize, "time",timestamp, "txt", sep="." )
      full_age_file = paste( full_age_dir,'/', age.file.name, sep='')
      
      write.csv( popAges, full_age_file, row.names=F)
        
}  
