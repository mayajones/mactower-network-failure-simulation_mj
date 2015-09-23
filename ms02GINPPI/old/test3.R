#20150623 blacklight

rm(list=ls())
debug = 0
mydir = "/brashear/hqin2/mactower-network-failure-simulation-master/ms02GINPPI"
setwd("/brashear/hqin2/mactower-network-failure-simulation-master/ms02GINPPI")
source('network.r')

net = read.table( "merged_PPIGIN_2014Jan20.tab", header=F, sep="\t", colClass = c("character", "character") )
head(net)
if(debug==9) { 
  #net = read.table('pair.tab',header=F) 
 net = net[1:90000,]
}

for( i in 400:500) {
 net.ms02 = ms02_singlerun( net, indebug=0 )
 cmnd = paste( "mkdir ", mydir, "/dipgin.ms02.output/", i, sep="")
 system( cmnd )
 outputname = paste( mydir, '/dipgin.ms02.output/', i, '/', "ms02_",i,".tab", sep="")
 write.table(net.ms02, outputname, quote=F, row.names=F, sep="\t")  #2014Feb 17
}

quit('no')



