#20150914 testing parallel jobs on greenfield
# R CMD BATCH ./example.R --args 2 35 red
rm(list=ls())

source("lifespan.r")
source("network.r")

#R -f file --args start end
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))
start = as.integer(args[1]); start; 
end = as.integer(args[2]); end; 
tag = args[3]
  
myhost = 'greenfield'  # 'byte' 'blacklight' 'mactower'
#myhost = 'byte'  # 'byte' 'blacklight' 'mactower'
#myhost = 'helen'  # 'byte' 'blacklight' 'helen'

mydir = "/crucible/mc48o9p/hqin2/mactower-network-failure-simulation-master/ms02GINPPI"
if (myhost == 'byte') {  mydir = "/Users/hqin/github/mactower-network-failure-simulation/ms02GINPPI"
} else if (myhost == 'helen') { mydir = "/Users/hqin/github/mactower-network-failure-simulation/ms02GINPPI";  
}

print(paste('tag is', tag, "  current dir is:", mydir))
#list.files()
getwd()

print(paste("Now setwd"))
setwd(mydir)
#list.files()

debug = 0; 
x = start:end 
outfile = paste( tag, start, end, "tab", sep='.')
write.csv( x, outfile, row.names=F)
      
    
