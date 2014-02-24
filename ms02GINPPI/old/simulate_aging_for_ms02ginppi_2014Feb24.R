#permuate merged yeast PPI+GIN

#2014 Feb 12, re-name function to ms02_singlerun
#2014 Jan 31, fixed a bug that inserted "NA" into new network. The bug seems to be caused by spliting the 
# arrays. I rewrote the spliting portion. 

#require(igraph)
rm(list=ls())
debug = 9
setwd("~/projects/0.ginppi.reliability.simulation/ms02GINPPI")
#set.seed(2014)

source("network.r")
require('flexsurv')
#source("/Users/hongqin/lib/R/lifespan.r")
source("lifespan.r")

#essential gene info
essenTb = read.csv('SummaryRegressionHetHom2013Oct29.csv', colClasses=rep('character', 9))




 
 
}


