source("lifespan.r")

#lifespan = round(rnorm(1000)*10, 0) + 30
#lifespan = lifespan[lifespan>1]

lifespan = round( rgompertz(0.001,0.2, 10000) + 0.5 )
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

plot( log10(tb$mortality.rate) ~ log10(tb$t), type='l'  ) #linear for Weibull
plot( log10(tb$mortality.rate) ~ tb$t, type='l') #linear for Gompertz

tb2 = tb2[-length(tb2[,1]), ]

summary(lm(log10(tb2$mortality.rate) ~ tb2$t ))

summary(lm(log10(tb2$mortality.rate) ~ log10(tb2$t) ))

require(flexsurv)
lifespanGomp = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
lifespanWeib = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'weibull')  

c(lifespanWeib$AIC, lifespanGomp$AIC, lifespanWeib$AIC - lifespanGomp$AIC )



