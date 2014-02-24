
setwd("~/Dropbox/Alexander_Ramp_Research/robustnessAA/current2012Fall/data")

###INPUT DATA INTO R


#Read lifespan table into R
lifespan = read.csv("lifespan.csv");

#Read growth fitness table into R
fit = read.csv("growth.fitness.hom.csv");
fit$orf = as.character( fit$orf ); 

#Read evolutionary distance table into R
Kdata = read.csv( "Sce.Spa.KaKs.csv");

#Input genetics interaction pairs
pairs = read.csv("pairs.csv");
gPairs = read.csv("sgadata_costanzo2009_stringentCutoff_101120.csv", header=F)
pairs = tmp[,c(1,3)]
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

#Read mutant table into R
scmd = read.table( "scmd.tab", sep="\t", header=T)

#Normalize the data by column
row.names(scmd) = as.character( scmd$name )
for( j in 2:502 ){
 scmd[,j] = ( scmd[,j] - mean(scmd[,j],na.rm=T) )/ sqrt( var( scmd[,j], na.rm=T))
};

#Calculate the standard deviation by row
for ( i in 1:4718){
scmd$stddev[i] = sqrt(var( t(scmd[i,2:502]), na.rm=T ) )}; 

#save.image("Matheson4_2012.2,00pm.RData")
#load("Matheson4_2012.2,00pm.RData")

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
lifespan$Ka = Sce.Smik$Ka[match(lifespan$ORF,data$orfname)];
m = lm(lifespan$RLS_Del_alpha ~ lifespan$Ka);
summary(m);

## GENETIC INTERACTIONS VS LIFELIFESPAN REGRESSION ANALYSIS ##
#calculate connecting degrees for proteins
ids = c(pairs$ORF1, pairs$ORF2);
degree = table( ids );

#Summarize the results in a dataframe 'net'
net = data.frame(degree);
str(net);
net$id = as.character( net$id);

#Match data and net
data$degree = net$Freq[match(data$orfname, net$id )]

#Match lifespan to Frequency
lifespan$Freq = net$Freq[match(lifespan$ORF, net$id )];

#Perform linear regression of lifespan vs evolutionary distance and summarize results
m= lm(lifespan$RLS_Del_alpha ~ lifespan$Freq + lifespan$YPE);
summary(m);

#Plot regression analysis
plot(lifespan$RLS_Del_alpha ~ lifespan$Freq, ylab="Replicative Lifespan", xlab="Genetic Interactions" );
abline( m, col="red");

###MORPHOLOGICAL PLASTICITY VS LIFESPAN REGRESSION ANALYSIS
#Match lifespan to mutant data
lifespan$stddev = scmd$stddev[match(lifespan$ORF, scmd$name)];

#Perform linear regression of lifespan to morphological plasticity and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$stddev);
summary(m);
plot(lifespan$RLS_Del_alpha ~ lifespan$stddev)

#Plot of regresssion analysis
plot(scmd$RLS_Del_alpha ~ scmd[,j]) ylab="Replicative Lifespan", xlab="Morphological Plasticity Robustness" );
abline( m, col="red");

###LIFESPAN MULTIPLE REGRESSION ANALYSIS
#Perform regression analysis and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$YPE + lifespan$Ka + lifespan$Freq + lifespan$stddev );
summary(m);

###ROBUST FACTOR MULTIPLE REGRESSION ANALYSIS
#Perform regression analysis and summarize results
m = lm(lifespan$YPE ~ lifespan$Ka + lifespan$Freq + lifespan$stddev);
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