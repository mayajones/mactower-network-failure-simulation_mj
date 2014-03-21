#plot mortality rate from simulated GINPPI aging.

rm(list=ls())
require(flexsurv)
setwd("~/github/mactower-network-failure-simulation/mathramp")
source("lifespan.r")

getwd()
list.files()

origAgeFiles = list.files(path="original_ginppit_failures")
origAgeFiles

list.files(path="ms02.gindip.failures/2/popages")
ms02AgeFiles = list.files(path='ms02.gindip.failures/1/popages', pattern='2000')
ms02AgeFiles

tb.ori = read.csv("original_ginppit_failures/cutoff.4.p.0.7.lambda.0.04.popsize.10000.time.2014Mar11_233122.txt")
tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.7.lambda.0.04.popsize.500.time.2014Feb27_153343.txt")

#tb.ori = read.csv("original_ginppit_failures/cutoff.4.p.0.9.lambda.0.04.popsize.10000.time.2014Mar16_091252.txt")
#tb.ms02 = read.csv("ms02.gindip.failures/1/popages/cutoff.4.p.0.9.lambda.0.04.popsize.500.time.2014Feb27_205447.txt")

summary(tb.ori)
summary(tb.ms02)
s.ori = calculate.s(tb.ori[,1])
s.ms02 = calculate.s(tb.ms02[,1])

ks.test( tb.ori[,1], tb.ms02[,1])
wilcox.test( tb.ori[,1], c(tb.ms02[,1],tb.ms02[,1], tb.ms02[,1], tb.ms02[,1]))

#First, plot of original network survival curve
plot( s.ori$s ~ s.ori$t, type='l',col='blue' )
lines( s.ms02$s ~ s.ms02$t, col='red')

plot( s.ori$s ~ s.ori$t, type='l',col='blue', log='x' )
lines( s.ms02$s ~ s.ms02$t, col='red')

#fitting Gompertz and Weibull
ori.gomp = flexsurvreg(formula = Surv(tb.ori[1:500,1]) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ms02.gomp = flexsurvreg(formula = Surv(tb.ms02[1:500,1]) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ori.weibull = flexsurvreg(formula = Surv(tb.ori[1:500,1]) ~ 1, dist = 'weibull') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
ms02.weibull = flexsurvreg(formula = Surv(tb.ms02[1:500,1]) ~ 1, dist = 'weibull') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
c(ori.gomp$AIC, ms02.gomp$AIC, ori.weibull$AIC, ms02.weibull$AIC)


lifespan = round( tb.ori[,1], digits=1)
table(lifespan)
tb = calculate.s(lifespan)
head(tb)
tb$ds=NA; tb$dt=NA

tb$dt[1] = tb$s[1]
tb$ds[1] = 1 - tb$s[1]
tb$mortality.rate[1] = tb$ds[1] / tb$dt[1]

for( j in 2:length(tb[,1])) {
  tb$ds[j] =  tb$s[j-1] - tb$s[j] 
  tb$dt[j] = -tb$t[j-1] + tb$t[j]
  tb$mortality.rate[j] = tb$ds[j] / ( tb$s[j] * tb$dt[j])
}
plot( tb$s ~ tb$t)
plot( tb$mortality.rate ~ tb$t, typ='l' )

plot( log10(tb$mortality.rate) ~ tb$t, type='p') #linear for Gompertz, semi-log plot

plot( log10(tb$mortality.rate) ~ log10(tb$t), type='p'  ) #linear for Weibull, log-log plot



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



