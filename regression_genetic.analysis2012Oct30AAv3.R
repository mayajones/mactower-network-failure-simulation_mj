
setwd("~/Dropbox/Alexander_Ramp_Research/robustnessAA/current2012Fall/data")

###INPUT DATA INTO R


#Read lifespan table into R
lifespan = read.csv("lifespan.csv");

#Read growth fitness table into R
fit = read.csv("growth.fitness.hom.csv");
fit$orf = as.character( fit$orf ); 

#Read evolutionary distance table into R
Kdata = read.csv( "Sce.Spa.KaKs.csv");

gPairs = read.csv("sgadata_costanzo2009_stringentCutoff_101120.csv", header=F)
names(gPairs) = c("ORF1", "Name1", "ORF2", "Name2", NA, NA, NA)
gPairs$ORF1 = as.character( gPairs$ORF1 )
gPairs$ORF2 = as.character( gPairs$ORF2 )


#Input protein interaction pairs
pairs = read.csv("pairs.csv");
colnames(pairs) = c("ORF1", "ORF2")
pairs$ORF1 = as.character( pairs$ORF1 );
pairs$ORF2 = as.character( pairs$ORF2 );

#Calculate connecting degrees for proteins
ids = c(pairs$ORF1, pairs$ORF2);
degree = table( ids );

#Summarize the results in a dataframe 'net'
net = data.frame(degree);
str(net);
net$id = as.character( net$id);

#Visual examination to double the merged results
head(data)
net[net$id=='YAL005C', ]
net[net$id=='YAL012W', ]

#Read cellular morpholy mutant table into R
scmd = read.table( "scmd.tab", sep="\t", header=T)

#Normalize the data by column
row.names(scmd) = as.character( scmd$name )
for( j in 2:502 ){
 scmd[,j] = ( scmd[,j] - mean(scmd[,j],na.rm=T) )/ sqrt( var( scmd[,j], na.rm=T))
};
head(scmd)

#Calculate the sigma, standard deviation by row
for ( i in 1:4718){
 scmd$stddev[i] = sqrt(var( t(scmd[i,2:502]), na.rm=T ) )
 scmd$mean[i] = mean( t(scmd[i,2:502]), na.rm=T ) 
}; 
scmd$CV = scmd$stddev / (scmd$mean - min(scmd$mean) + 0.01) #make sure CV is positive
head(scmd[,c("CV", "stddev", "mean")])
summary(scmd[,c("CV", "stddev", "mean")])

save.image("../output/AA2012Oct30-V1.RData")
load("../output/AA2012Oct30-V1.RData")

###FITNESS VS LIFESPAN REGRESSION ANALYSIS

#Determine which fitness growth medium has the best relationship to RLS
lifespan$YPD = fit$YPD[match(lifespan$ORF, fit$orf)];
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$YPD));

lifespan$YPDGE = fit$YPDGE[match(lifespan$ORF, fit$orf)];
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$YPDGE));

lifespan$YPG = fit$YPG[match(lifespan$ORF, fit$orf)];
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$YPG));

lifespan$YPE = fit$YPE[match(lifespan$ORF, fit$orf)];
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$YPE));

lifespan$YPL = fit$YPL[match(lifespan$ORF, fit$orf)];
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$YPL));

##The best fitness vs lifespan data was for YPE because the R squared value was the strongest and the p value was smallest.

#Plot YPE fitness vs lifespan
plot(lifespan$RLS_Del_alpha ~ lifespan$YPE, ylab="Replicative Lifespan", xlab="Fitness Robustness");
abline( m, col="red");


###EVOLUTIONARY DISTANCE VS LIFESPAN REGRESSION ANALYSIS
#Match lifespan to evolutionary distance
lifespan$Ka = data$Ka[match(lifespan$ORF,data$orfname)];

#Perform linear regression of lifespan vs evolutionary distance and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$Ka);
summary(m);

#Plot regression analysis
plot(lifespan$RLS_Del_alpha ~ lifespan$Ka, ylab="Replicative Lifespan", xlab="Evolutionary Distance(K)");
abline( m, col="red");

###NEW CODE###
#SCE VS.SMIK = Evolutionary distance ##
Sce.Smik = read.csv("Sce.Smik.csv");
lifespan$Ka = Sce.Smik$Ka[match(lifespan$ORF,data$orfname)];
m = lm(lifespan$RLS_Del_alpha ~ lifespan$Ka);
summary(m);
#SCE VS. SBAY = Evolutionary distance ##
Sce.Sbay = read.csv("Sce.Sbay.csv");
lifespan$Ka = Sce.Sbay$Ka[match(lifespan$ORF,data$orfname)];
m = lm(lifespan$RLS_Del_alpha ~ lifespan$Ka);
summary(m);

## GENETIC INTERACTIONS VS LIFELIFESPAN REGRESSION ANALYSIS ##
#calculate connecting degrees for proteins
gids = c(gPairs$ORF1, gPairs$ORF2);
gdegree = table( gids );

#Summarize the results in a dataframe 'net'
gnet = data.frame(gdegree);
str(gnet);
gnet$gids = as.character( gnet$gids);


#Match Kdata and gnet
Kdata$gDegree = gnet$Freq[match(Kdata$orfname, gnet$gids )]
head(Kdata)

#Match lifespan to Frequency
lifespan$gDegree = gnet$Freq[match(lifespan$ORF, gnet$gids )];
lifespan$pDegree =  net$Freq[match(lifespan$ORF,  net$ids )];
head(lifespan)

#Perform linear regression of lifespan vs evolutionary distance and summarize results
summary( lm(lifespan$RLS_Del_alpha ~                 lifespan$pDegree + lifespan$YPE))
summary( lm(lifespan$RLS_Del_alpha ~ lifespan$gDegree                 + lifespan$YPE))
summary( lm(lifespan$RLS_Del_alpha ~ lifespan$gDegree + lifespan$pDegree + lifespan$YPE))
summary( lm(lifespan$gDegree ~  lifespan$YPE))
summary( lm(lifespan$pDegree ~  lifespan$YPE))


#Plot regression analysis
plot(lifespan$RLS_Del_alpha ~ lifespan$gDegree, ylab="Replicative Lifespan", xlab="Genetic Interactions" );
abline( m, col="red");

# add CV for fitness plasticity
head(fit)
fitsub = fit[, c(9:17)]
head(fitsub)

fit$fitsd = apply(fitsub, 1, FUN= sd)
fit$fitmean = apply(fitsub, 1, FUN=mean)
fit$CV = fit$fitsd / fit$fitmean

lifespan$fit.CV<- fit$CV[match(lifespan$ORF, fit$orf)]
summary(lm( 1/ lifespan$fit.CV ~ lifespan$RLS_Del_alpha))
summary(lm( sqrt(1/ lifespan$fit.CV) ~ lifespan$RLS_Del_alpha))

#lifespan$sqrt1.over.fitCV = sqrt(1/ lifespan$fit.CV)
summary(lm( lifespan$RLS_Del_alpha ~  sqrt(1/ lifespan$fit.CV) + lifespan$pDegree)) #****** which reference???? Biological Implications of the Weibull and Gompertz Models of Aging

save.image("AA2012Oct30.RData")
load("AA2012Oct30.RData")

###MORPHOLOGICAL PLASTICITY VS LIFESPAN REGRESSION ANALYSIS
#Match lifespan to mutant data
lifespan$scmdstddev = scmd$stddev[match(lifespan$ORF, scmd$name)];
lifespan$scmdCV = scmd$CV[match(lifespan$ORF, scmd$name)];
lifespan$scmdMean = scmd$mean[match(lifespan$ORF, scmd$name)];

#Perform linear regression of lifespan to morphological plasticity and summarize results
summary(lm(lifespan$RLS_Del_alpha ~ lifespan$scmdstddev))
plot(lifespan$RLS_Del_alpha ~ lifespan$scmdstddev)

summary(lm(lifespan$RLS_Del_alpha ~ sqrt(1/lifespan$scmdstddev)))
summary(lm(lifespan$RLS_Del_alpha ~ sqrt(lifespan$scmdstddev))) #extreem small p

summary(lm(lifespan$RLS_Del_alpha ~ lifespan$scmdMean)) #p=0.0019

summary(lm(lifespan$scmdstddev ~ lifespan$scmdMean)) ## highly correlated mean and stddev in morphology??

summary(lm(lifespan$RLS_Del_alpha ~ lifespan$scmdCV))  #does mean and stddev offset each other? 
#summary(lm(lifespan$RLS_Del_alpha ~ sqrt(lifespan$scmdCV))


#Plot of regresssion analysis
plot(scmd$RLS_Del_alpha ~ scmd[,j], ylab="Replicative Lifespan", xlab="Morphological Plasticity Robustness" );
abline( m, col="red");

###LIFESPAN MULTIPLE REGRESSION ANALYSIS
#Perform regression analysis and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$YPE + lifespan$Ka + lifespan$Freq + lifespan$stddev );
summary(m);

###ROBUST FACTOR MULTIPLE REGRESSION ANALYSIS
#Perform regression analysis and summarize results
#m = lm(lifespan$YPE ~ lifespan$Ka + lifespan$Freq + lifespan$stddev); #!!!!!!!!!!!!Error lifespan$Freq
m = lm(lifespan$YPE ~ lifespan$Ka + lifespan$pDegree + lifespan$stddev);
summary(m);

m = lm(lifespan$YPE ~ lifespan$Ka + lifespan$gDegree + lifespan$stddev); 
summary(m);

###FITNESS, MORPHOLOGICAL PLASTICITY & LIFESPAN MULTIPLE REGRESSION
#Perform regression analysis and summarize results when RLS is held constant
m = lm(lifespan$RLS_Del_alpha ~ lifespan$stddev + lifespan$YPE);
summary(m);

#Perform regression analysis and summarize results when morphological plasticity is held constant
m = lm(lifespan$stddev ~ lifespan$YPE + lifespan$RLS_Del_alpha);
summary(m);

#Perform regression analysis and summarize results when fitness is held constant
m = lm(lifespan$YPE ~ lifespan$stddev + lifespan$RLS_Del_alpha);
summary(m);