#permuate merged yeast PPI+GIN

#2015 June 25. Commandline pass parameters
#2014 April 9. Rerun using 2014 April 9 version 'network.r' with correct ms02 function.
#2014 Feb 12, re-name function to ms02_singlerun
#2014 Jan 31, fixed a bug that inserted "NA" into new network. The bug seems to be caused by spliting the 
# arrays. I rewrote the spliting portion. 

rm(list=ls())

#R -f file --args start end
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))
start = as.integer(args[1])
end = as.integer(args[2])

debug = 0
myhost = 'greenfield'  # 'byte' 'blacklight' 'mactower'
#myhost = 'byte'  # 'byte' 'blacklight' 'mactower'

mydir = "/crucible/mc48o9p/hqin2/mactower-network-failure-simulation-master/ms02GINPPI"
if (myhost == 'byte') {
  mydir = "~/github/mactower-network-failure-simulation/ms02GINPPI"
}  

setwd(mydir)
list.files()
#set.seed(2014)
source('network.r')

#net = read.table("repeat.tab")
#write.table(pairs, "merged_PPIGIN_2014Jan20.tab", quote=F, row.names=F, col.names=F, sep='\t')
net = read.table( "merged_PPIGIN_2014Jan20.tab", header=F, sep="\t", colClass = c("character", "character") )
head(net)
if(debug==9) { 
  #net = read.table('pair.tab',header=F) 
 net = net[1:90000,]
}

for( i in start:end) {
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





