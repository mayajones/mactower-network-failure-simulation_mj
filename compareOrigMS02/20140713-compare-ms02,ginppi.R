#compare agine dynamcis between original and ms02 GINDIP. 
# plot mortality rate ~time from simulated GINPPI aging.
# binomial aging model fitting

rm(list=ls())
setwd("~/github/mactower-network-failure-simulation/compareOrigMS02")
require(flexsurv)
getwd()
list.files()

source("lifespan.20140711.R")

ms02path = 'ms02.gindip.failures'
list.files(, path=ms02path)

ms02fullnames = NA;
for(i in 1:20){
  mypath = paste( ms02path, '/', i, '/popages', sep='' )
  myfiles = list.files(, path=mypath)
  myfullnames = paste( mypath, '/', myfiles, sep='')
  ms02fullnames = c(ms02fullnames, myfullnames)  
}
ms02fullnames = ms02fullnames[-1]

origAgeFiles = list.files(path="ori_ginppit_ages")
origAgeFiles

mypattern = "p.0.8.lambda.0.01.popsiz"
ms02files= ms02fullnames[grep(mypattern, ms02fullnames)]
ms02output = data.frame(ms02files); 
ms02output$mean = NA; ms02output$stddev = NA; ms02output$CV=NA; 
ms02output$R=NA; ms02output$t0=NA; ms02output$n=NA;
ms02output$Rflex = NA; ms02output$Gflex = NA; 

for( i in 1:length(ms02files)) {
  tb = read.csv(ms02files[i])
  summary(tb); names(tb) = c('rls'); tb$rls = round(tb$rls+0.5, digits=0)
  ms02output$mean[i]=mean(tb$rls)
  ms02output$stddev[i]= sd(tb$rls)
  ms02output$CV[i] = sd(tb$rls) / mean(tb$rls)
  
  fitGomp = flexsurvreg( formula=Surv(tb$rls) ~ 1, dist='gompertz')
  Rhat = fitGomp$res[2,1]; Ghat = fitGomp$res[1,1] 
  ms02output$Rflex[i] = Rhat; 
  ms02output$Gflex[i]=Ghat; 
  
  #G = (n-1) / t0, so t0 = (n-1)/G
  nhat = 4; t0= (nhat-1)/Ghat; t0
  fitBinom = optim ( c(Rhat, t0, nhat),  llh.binomialMortality.single.run , lifespan=tb$rls, 
                     method="L-BFGS-B", lower=c(1E-10,0.05,1E-1), upper=c(10,100,50) );  
  
  ms02output[i, c("R", "t0", "n")] = fitBinom$par
  
  tb2 = calculate.mortality.rate(tb$rls)
  plot( tb2$s ~ tb2$t, main=ms02files[i])
  plot( tb2$mortality.rate ~ tb2$t , main=ms02files[i])
  plot( tb2$mortality.rate ~ tb2$t, log='y', main= paste( 'linear for Gompertz', main=ms02files[i]) )
  plot( tb2$mortality.rate ~ tb2$t, log='yx', main=paste( 'linear for Weibull',main=ms02files[i]) )
}#i file loop

ms02output

oripath = "ori_ginppit_ages"
#origAgeFiles = list.files(path= oripath)
origAgeFiles
origFile = origAgeFiles[grep(mypattern, origAgeFiles)]
tb.ori = read.csv(paste( oripath, '/',origFile[2],sep='')  )
tb.ori$rls = round( tb.ori[,1]+0.5, digit=0)

oriGomp = flexsurvreg( formula=Surv(tb.ori$rls) ~ 1, dist='gompertz')
Rori = oriGomp$res[2,1]; Gori = oriGomp$res[1,1] 
nori = 4; t0= (nori-1)/Gori; t0
oriBinom = optim ( c(Rori, t0, nori),  llh.binomialMortality.single.run , lifespan=tb.ori$rls, 
                   method="L-BFGS-B", lower=c(1E-10,0.05,1E-1), upper=c(10,100,50) );  

cbind( mean(tb.ori$rls), median(ms02output$mean))
cbind( sd(tb.ori$rls), median(ms02output$stddev))
cbind( sd(tb.ori$rls)/mean(tb.ori$rls), median(ms02output$CV))
cbind( Rori, median(ms02output$Rflex))
cbind( oriBinom$par[1], median(ms02output$R))
cbind( oriBinom$par[2], median(ms02output$t0))
cbind( oriBinom$par[3], median(ms02output$n)) #ms02 does change n
cbind( Gori, median(ms02output$Gflex))

cv.ori =  sd(tb.ori$rls)/mean(tb.ori$rls)
hist(ms02output$CV, br=15)
arrows(cv.ori, 2.5, cv.ori, 1.1, lwd=2, col='red', xlab="t0" )

hist(ms02output$Rflex)
hist(ms02output$t0)
# arrows(oriBinom$par[2], 5, oriBinom$par[2], 1, lwd=2, col='red', xlab="t0" )

s.ori = calculate.mortality.rate(tb.ori$rls)
plot( s.ori$mortality.rate ~ s.ori$t, type='l', col='red', log='y', xlim=c(40,70), lwd=3)

for( i in 1:10) {
tb = read.csv(ms02files[i]); tb$rls = round(tb[,1]+0.5, digits=0)
s2 = calculate.mortality.rate(tb$rls )
lines(s2$mortality.rate ~ s2$t, col='gray', lwd=1)
}

lines(s2$s ~ s2$t)

plot( s.ori$mortality.rate ~ s.ori$t, type='l', col='red')
lines(s2$mortality.rate ~ s2$t, col='blue')

plot( s.ori$mortality.rate ~ s.ori$t, type='l', col='red', log='y')
lines(s2$mortality.rate ~ s2$t, col='blue')


quit('no')
##########END 20140713










#ks.test( tb.ori[,1], tb.ms02[,1])
t.test( tb.ori[,1], tb.ms02[,1])
d1 = density(tb.ori[,1])
d2 = density(tb.ms02[,1])
plot( d1$y ~ d1$x, type='l', col='blue', xlab='age', ylab='density', lwd=4, ylim=c(0,0.035), xlim=c(0,120))
lines( d2$y ~ d2$x, col='red', lwd=4)
legend(80,0.03, c("original","random"), lwd=c(4,4), col=c('blue','red'))

plot( d1$y ~ d1$x, type='l', col='blue', xlab='age', ylab='density', lwd=4)
lines( d2$y ~ d2$x, col='red', lwd=4)

#First, plot of original network survival curve
plot( s.ori$s ~ s.ori$t, type='l',col='blue' )
lines( s.ms02$s ~ s.ms02$t, col='red')

plot( s.ori$s ~ s.ori$t, type='l',col='blue', log='x' )
lines( s.ms02$s ~ s.ms02$t, col='red')

#fitting Gompertz and Weibull
ori.gomp = flexsurvreg(formula = Surv(tb.ori[1:2000,1]) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ms02.gomp = flexsurvreg(formula = Surv(tb.ms02[1:2000,1]) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ori.weibull = flexsurvreg(formula = Surv(tb.ori[1:2000,1]) ~ 1, dist = 'weibull') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ms02.weibull = flexsurvreg(formula = Surv(tb.ms02[1:2000,1]) ~ 1, dist = 'weibull') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
c(ori.gomp$AIC, ori.weibull$AIC, ms02.gomp$AIC, ms02.weibull$AIC)
c((ori.gomp$AIC - ori.weibull$AIC), (ms02.gomp$AIC - ms02.weibull$AIC) )


#calculate mortality rate
calculate.mortality.rate = function( lifespan ){
  tb = calculate.s(lifespan)
  tb$ds=NA; tb$dt=NA
  #first point
  tb$dt[1] = tb$t[2]
  tb$ds[1] = 1 - tb$s[1]
  tb$mortality.rate[1] = tb$ds[1] / tb$dt[1]
  
  for( j in 2:length(tb[,1])) {
    tb$ds[j] =  tb$s[j-1] - tb$s[j] 
    tb$dt[j] = -tb$t[j-1] + tb$t[j]
    tb$mortality.rate[j] = tb$ds[j] / ( tb$s[j] * tb$dt[j])
  }
  return(tb)
} #end of calculate.mortality.rate()

tb.ori[1:2000,1])
tb.ms02[1:2000,1])

out.ori = calculate.mortality.rate(round( tb.ori[,1], digits=1))
out.ms02 = calculate.mortality.rate(round( tb.ms02[,1], digits=1))

plot( out.ori$mortality.rate ~ out.ori$t, type='l', log='y'  ) #linear for Gompertz, semi log  plot
plot( out.ori$mortality.rate ~ out.ori$t, type='l', log='xy'  ) #linear for Weibull, log-log plot

summary( lm( log(out.ori$mortality.rate)~out.ori$t ), na.rm=T )

plot( out.ms02$mortality.rate ~ out.ms02$t, type='l', log='y'  ) #linear for Gompertz, semi log  plot
plot( out.ms02$mortality.rate ~ out.ms02$t, type='l', log='xy'  ) #linear for Weibull, log-log plot



#########end

lifespan = round( tb.ori[,1], digits=1)
table(lifespan)
tb = calculate.s(lifespan)
head(tb)
tb$ds=NA; tb$dt=NA
tb$dt[1] = tb$t[2] #20130321 correct a bug tb$s -> tb$t
tb$ds[1] = 1 - tb$s[1]
tb$mortality.rate[1] = tb$ds[1] / tb$dt[1]

for( j in 2:length(tb[,1])) {
  tb$ds[j] =  tb$s[j-1] - tb$s[j] 
  tb$dt[j] = -tb$t[j-1] + tb$t[j]
  tb$mortality.rate[j] = tb$ds[j] / ( tb$s[j] * tb$dt[j])
}
plot( tb$s ~ tb$t)
plot( tb$mortality.rate ~ tb$t, typ='l', log='y' )

plot( log10(tb$mortality.rate) ~ tb$t, type='l') #linear for Gompertz, semi-log plot
plot( log10(tb$mortality.rate) ~ log10(tb$t), type='l'  ) #linear for Weibull, log-log plot

#the first point was very low because it is a bug.
#and last point of mortality rate are very low, presumbaly due to boundary effect.?!



###############old file

title(datafiles[i])
title(datafiles[i])

tb2 = tb 
tb2 = tb2[-length(tb2[,1]), ]


summary(lm(log10(tb2$mortality.rate) ~ tb2$t ))

summary(lm(log10(tb2$mortality.rate) ~ log10(tb2$t) ))

require(flexsurv)
lifespanGomp = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
lifespanWeib = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'weibull')  

c(lifespanWeib$AIC, lifespanGomp$AIC, lifespanWeib$AIC - lifespanGomp$AIC )



