

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
