# generate gompertz random numbers
# fit with flexsurv Gompertz and weibull models

#inverse of gompertz CDF
# see http://hongqinlab.blogspot.com/2013/06/median-lifespan-of-2-parameter-gompertz.html
inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }

#see generate random number with a given distribution
# http://hongqinlab.blogspot.com/2013/12/generate-gompertz-random-numbers.html
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
  x.gompertz = inverse.gomp.CDF(0.001,0.2, x.uniform)
  return(x.gompertz)
}
rgompertz(0.001,0.2,100)

x.uniform = runif(500)
hist(x.uniform)

x.gompertz = inverse.gomp.CDF(0.001,0.2, x.uniform)
hist(x.gompertz)
summary(x.gompertz)

source("lifespan.r")
tb = calculate.s(x.gompertz)
plot(tb$s ~ tb$t)

require(flexsurv)
lifespan = x.gompertz
lifespanGomp = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'gompertz')  
lifespanWeib = flexsurvreg(formula = Surv(lifespan) ~ 1, dist = 'weibull')  

c(lifespanWeib$AIC, lifespanGomp$AIC, lifespanWeib$AIC - lifespanGomp$AIC )

