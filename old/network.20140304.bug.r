

single_network_failure = function(lambda, p, pairs, runningORFs) {
  # single network failure simulation
  # lambda: exponential constant failure rate for edges
  # pairs: network in pairwide format
  # runningORFs: GooddEssentialORFsPPI  
  
  inpairs = pairs[,1:2] #bookkeeping  
  names(inpairs) = c('id1','id2')
  
  #stochasticity into pairs   
  inpairs$active = runif(length(inpairs[,1]))  #uniform
  # tmp = pairs$active > 1-p
  # table(tmp) / length(tmp)  ; #double-check, very good. 
  
  inpairs$age = rexp( length(inpairs[,1]), rate=lambda )  #exponential ages for pairs
  inpairs$age = ifelse(inpairs$active > (1-p), inpairs$age, NA ) #if not active, intxn is excluded. 
  #pairs$age = ifelse(pairs$active > (1-p), pairs$age, 0 )  # in what situations, can non-ative intxn be treat as 0-age?
  
  ModuleTb = data.frame(runningORFs) #buffer for module ages    
  #loop every essential genes to identify the module age
  for (i in 1:length(runningORFs)) {
    myORF = runningORFs[i]
    pos1 = grep(myORF, inpairs$id1)
    pos2 = grep(myORF, inpairs$id2)  #id1,2 to ORF1,2 is a really bad choice. 
    if( length( c(pos1,pos2))>=1 ) {
      ModuleTb$age.m[i] = max( inpairs$age[c(pos1,pos2)], na.rm=T )   #maximal intxn age -> module age
    } else {
      ModuleTb$age.m[i] = NA; 
    }
  }
  #head(ModuleTb); 
  summary(ModuleTb)
  ModuleTb$age.m[ ModuleTb$age.m== -Inf] = 0; #dead births occur when links are not active
  currentNetworkAge = min(ModuleTb$age.m)
}



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
