#2014 April 9. Rerun using 2014 April 9 version 'network.r' with correct ms02 function.

#permuate merged yeast PPI+GIN

#2014 Feb 12, re-name function to ms02_singlerun
#2014 Jan 31, fixed a bug that inserted "NA" into new network. The bug seems to be caused by spliting the 
# arrays. I rewrote the spliting portion. 



#require(igraph)
rm(list=ls())
debug = 0
setwd("~/projects/0.ginppi.reliability.simulation/ms02GINPPI")
#set.seed(2014)

#permute.pairs.wo.selfpairs = function( inpairs,  ncycles=10, debug=1 ) {  
ms02_singlerun = function( inpairs,  ncycles=10, indebug=0 ) { # Renamed, 2014 Feb 12
  if (ncycles >= 1 ) {
    if(indebug>0) {
      print(paste('ncycles=', ncycles))
    }
    longids = c(as.character(inpairs[,1]), as.character(inpairs[,2]) )
    longids = sample(longids)
    len = length(inpairs[,1])
    newpairs = data.frame( cbind( longids[1:len], longids[(len+1): (2*len)]) )
    names(newpairs) = c('id1', 'id2')
    newpairs$id1 = as.character( newpairs$id1)
    newpairs$id2 = as.character( newpairs$id2)    
    newpairs$selfpairs = ifelse( newpairs$id1 == newpairs$id2, 1, 0 )
    self.tb = newpairs[ newpairs$selfpairs==1, ]
    nonself.tb = newpairs[newpairs$selfpairs==0, ]
    if(indebug>0) {
      print(paste("===selfpairs===="),NULL)
      print(self.tb)
      print(paste("================="),NULL)
    }
    if( length(self.tb[,1])>=1 ) {
      if ( ncycles == 0) { 
        #return (c(NA,NA, NA) );
        print(paste("ncycles reached zero, ncycles"),ncycles)
        print(paste("Abort!"),NULL)
        stop; 
      } else {
        ncycles = ncycles - 1
        splitPos = round( length(self.tb[,1]) * sqrt(ncycles) ) + 5  #2014Jan31 change 
        splitPos = min( splitPos, (length(nonself.tb[,1])-1 ) )
        selectedpairs = rbind(self.tb,  nonself.tb[1: splitPos, ] )
        restpairs = nonself.tb[ (splitPos + 1): length(nonself.tb[,1]), ]
        #return( rbind(restpairs, permute.pairs.wo.selfpairs(selectedpairs, ncycles)))
        return( rbind(restpairs, ms02_singlerun(selectedpairs, ncycles)))  #2014 Feb 12
      }
    } else {  
      return (newpairs)
    }
  } else {
    return( c(NA,NA,NA )) 
  }
}

#net = read.table("repeat.tab")
#write.table(pairs, "merged_PPIGIN_2014Jan20.tab", quote=F, row.names=F, col.names=F, sep='\t')
net = read.table( "merged_PPIGIN_2014Jan20.tab", header=F, sep="\t", colClass = c("character", "character") )
head(net)
if(debug==9) { 
  #net = read.table('pair.tab',header=F) 
 net = net[1:90000,]
}

for( i in 1:100) {
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





