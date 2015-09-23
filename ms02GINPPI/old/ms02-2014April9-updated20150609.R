#20150623 blacklight
#2014 April 9. Rerun using 2014 April 9 version 'network.r' with correct ms02 function.

#permuate merged yeast PPI+GIN

#2014 Feb 12, re-name function to ms02_singlerun
#2014 Jan 31, fixed a bug that inserted "NA" into new network. The bug seems to be caused by spliting the 
# arrays. I rewrote the spliting portion. 


rm(list=ls())
debug = 0
setwd("/brashear/hqin2/mactower-network-failure-simulation/ms02GINPPI")
source('network.r')

#net = read.table("repeat.tab")
#write.table(pairs, "merged_PPIGIN_2014Jan20.tab", quote=F, row.names=F, col.names=F, sep='\t')
net = read.table( "merged_PPIGIN_2014Jan20.tab", header=F, sep="\t", colClass = c("character", "character") )
head(net)
if(debug==9) { 
  #net = read.table('pair.tab',header=F) 
 net = net[1:90000,]
}

for( i in 100:105) {
 net.ms02 = ms02_singlerun( net, indebug=0 )
 cmnd = paste( "mkdir dipgin.ms02.output/", i, sep="")
 system( cmnd )
 outputname = paste( 'dipgin.ms02.output/', i, '/', "ms02_",i,".tab", sep="")
 write.table(net.ms02, outputname, quote=F, row.names=F, sep="\t")  #2014Feb 17
}

#do they have the same degree?
#t1 = table(c(net[,1],net[,2]))
#t2 = table(c(net.ms02[,1],net.ms02[,2]))
#comp <- t1 == t2
#table(comp)
#tf = comp[comp==F]; tf
#t1[names(tf)[1]]
#t1[names(tf)]
#t2[names(tf)]





