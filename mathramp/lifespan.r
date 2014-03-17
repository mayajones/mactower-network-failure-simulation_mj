# source("/Users/hongqin/lib/R/lifespan.r")


##########################################
# single network failure simulation
# lambda: exponential constant failure rate for edges
# pairs: network in pairwide format
# runningORFs: GooddEssentialORFsPPI  

single_network_failure = function(lambda, p, pairs, runningORFs) {  
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


##########################################
#inverse of gompertz CDF
# see http://hongqinlab.blogspot.com/2013/06/median-lifespan-of-2-parameter-gompertz.html
#inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }

#see generate random number with a given distribution
# http://hongqinlab.blogspot.com/2013/12/generate-gompertz-random-numbers.html

#generate Gompertz random numbers
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
  x.gompertz = inverse.gomp.CDF(0.001,0.2, x.uniform)
  return(x.gompertz)
}
rgompertz(0.001,0.2,100)

##########################################
####
 #calculate s
 calculate.s = function( lifespan ){
 	 my.data = sort( lifespan[!is.na(lifespan)] );
   tmptb = table( my.data )
    for( i in 2:length(tmptb)) {
    	tmptb[i] = tmptb[i-1] + tmptb[i]    		} 	 
 	  tot = length(my.data)
 	 tmptb = tmptb / tot; 
 	 s = 1 - tmptb
 	 #list( s=s, t=unique(my.data));
 	 ret = data.frame( cbind(s, unique(my.data)));
 	 names(ret) = c("s", "t");
 	 ret;
 	}

##########################################
##### log likelihood function, simple gompertz mortality model
   #s = exp( (I/G) *(1 - exp(G* my.data)) )  ;
   #m = I * exp( G * my.data ) ;   
 llh.G.single.run <- function( IG, lifespan ) {
   I = IG[1]; G = IG[2];  
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data)) 
   #log_m = log( I * exp(G*my.data) + M); 
   log_m = log( I) + G*my.data; 
   my.lh = sum(log_s)  + sum(log_m);
   print (c(IG, my.lh) ); #trace the convergence
   ret = - my.lh # because optim seems to minimize
 }

##### log likelihood function, Gompertz-Makeham mortality model 
   #s = exp( (I/G) *(1 - exp(G* my.data)) - M*t )  ;
   #m = I * exp( G * my.data ) + M;   
 llh.GM.single.run <- function( IGM, lifespan ) {
   I = IGM[1]; G = IGM[2]; M = IGM[3]; 
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data)) - M*my.data
   log_m = log( I * exp(G*my.data) + M); 
   #log_m = log( I) + G*my.data; 
   my.lh = sum(log_s)  + sum(log_m);
   print (c(IGM, my.lh) ); #trace the convergence
   ret = - my.lh # because optim seems to minimize
 }

##########################################
   #s = exp( (I/G) *(1 - exp(G* my.data)) - M*t )  ;
 GM.s = function( IGM, t ) {
 	 I = IGM[1]; G = IGM[2]; M = IGM[3]; 
   log_s = (I/G) *(1 - exp(G* t)) - M*t
   ret  = exp( log_s )
 	}

#####
 G.s = function( IGM, t ) {
   I = IGM[1]; G = IGM[2]; M = 0; 
   log_s = (I/G) *(1 - exp(G* t)) - M*t
   ret  = exp( log_s )
 }
 
 #####
 
