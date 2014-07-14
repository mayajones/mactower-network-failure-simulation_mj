#compare agine dynamcis between original and ms02 GINDIP. 
# plot mortality rate ~time from simulated GINPPI aging.
# Gompertz and Weibull model fitting

rm(list=ls())
require(flexsurv)
setwd("~/github/mactower-network-failure-simulation/mathramp")
source("lifespan.r")

getwd()
list.files()

origAgeFiles = list.files(path="original_ginppit_failures2")
origAgeFiles

#list.files(path="ms02.gindip.failures/1/popages", pattern='2000')
ms02AgeFiles = list.files(path='ms02.gindip.failures/1/popages', pattern='2000')
#ms02AgeFiles = list.files(path='ms02-sim')
ms02AgeFiles

#tb.ori = read.csv("original_ginppit_failures/cutoff.4.p.0.7.lambda.0.04.popsize.10000.time.2014Mar11_233122.txt")
#tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.7.lambda.0.04.popsize.500.time.2014Feb27_153343.txt")

#tb.ori = read.csv("original_ginppit_failures2/cutoff.4.p.0.6.lambda.0.01.popsize.2000.time.2014Mar20_114812.txt")
#tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.6.lambda.0.01.popsize.2000.time.2014Mar18_190327.txt")

tb.ori = read.csv("original_ginppit_failures2/cutoff.4.p.0.7.lambda.0.01.popsize.2000.time.2014Mar21_000632.txt")
tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.7.lambda.0.01.popsize.2000.time.2014Mar19_074024.txt")

#tb.ori = read.csv("original_ginppit_failures2/cutoff.4.p.0.6.lambda.0.01.popsize.2000.time.2014Mar20_114812.txt")
#tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.6.lambda.0.01.popsize.2000.time.2014Mar18_190327.txt")

#for perfect networks
#tb.ori = read.csv("original_ginppit_failures2/cutoff.4.p.1.lambda.0.01.popsize.2000.time.2014Mar22_075849.txt")
#tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.1.lambda.0.01.popsize.2000.time.2014Mar20_155126.txt")
#tb.ms02= read.csv("ms02-sim/ms02_p1_lambda.0.01_pop4500.csv", header=F)

summary(tb.ori)
summary(tb.ms02)
s.ori = calculate.s(tb.ori[,1])
s.ms02 = calculate.s(tb.ms02[,1])

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



