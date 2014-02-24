###INPUT DATA INTO R
#Read lifespan table into R
lifespan = read.csv("lifespan.csv");

## will CLS data also be helpful? 

#Read growth fitness table into R
# http://www-deletion.stanford.edu/YDPM/Download_Data/Yeast_Deletion_Project/Regression_Tc1_hom.txt
fit = read.csv("growth.fitness.hom.csv");
fit$orf = as.character( fit$orf ); 
#There are new and alternative fitness data available

#Read evolutionary distance table into R
data = read.csv( "Sce.Spa.KaKs.csv");
# alternative KaKs should be used for cross-validation

#Input protein interaction pairs
pairs = read.csv("pairs.csv");
pairs$ORF1 = as.character( pairs$ORF1 );
pairs$ORF2 = as.character( pairs$ORF2 );
# we also use genetic interaction data (aka synthetic-lethal data)

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

# google SCMD and yeast
#Read morphology mutant data into R
scmd = read.table( "scmd.tab", sep="\t", header=T)

#Normalize the data by column by Z-score
row.names(scmd) = as.character( scmd$name )
for( j in 2:502 ){
 scmd[,j] = ( scmd[,j] - mean(scmd[,j],na.rm=T) )/ sqrt( var( scmd[,j], na.rm=T))
};

#Calculate the standard deviation by row (a proxy for morphology plasticity)
for ( i in 1:4718){
 scmd$stddev[i] = sqrt(var( t(scmd[i,2:502]), na.rm=T ) )
}; 

save.image("Metheson2012March27.2,40pm.RData")
load("Metheson2012March27.2,40pm.RData")

# Please go find, edit, and load protein half-life data into R. 2006 Belle...?
half.life = read.csv("protein half life.csv");

# Please see the referecen in the Dhami paper
# ... ... 

# Bonus point, mRNA half-life data


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
pdf("_RLS-YPE.pdf", width=5, height=5)
plot(lifespan$RLS_Del_alpha ~ lifespan$YPE, ylab="Replicative Lifespan", xlab="Fitness Robustness");
m = lm(lifespan$RLS_Del_alpha ~ lifespan$YPE )
abline( m, col="red");
summary(m)
text(0.5, 35, "Rsq= 0.0632")
text(0.5, 33, "p=0.0001084")
dev.off()
#YPE is ethenol and YPD is gluose. 
# Etnenol and gluose -> mitochondria metabolic state -> H2O2, O2*, etc. 

###EVOLUTIONARY DISTANCE VS LIFESPAN REGRESSION ANALYSIS
#Match lifespan to evolutionary distance
lifespan$Ka = data$Ka[match(lifespan$ORF,data$orfname)];

#Perform linear regression of lifespan vs evolutionary distance and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$Ka);
summary(m);

#Plot regression analysis
plot(lifespan$RLS_Del_alpha ~ lifespan$Ka, ylab="Replicative Lifespan", xlab="Evolutionary Distance(K)");
abline( m, col="red");


###PROTEIN INTERACTIONS VS LIFESPAN REGRESSION ANALYSIS
#calculate connecting degrees for proteins
ids = c(pairs$ORF1, pairs$ORF2);
degree = table( ids );

#Summarize the results in a dataframe 'net'
net = data.frame(degree);
str(net);
net$id = as.character( net$id);

#Match data and net
data$degree = net$Freq[match(data$orfname, net$id )]

#Match lifespan to Frequency (number of interaction per protein)
lifespan$Freq = net$Freq[match(lifespan$ORF, net$id )];

#Perform linear regression of lifespan vs evolutionary distance and summarize results
m= lm(lifespan$RLS_Del_alpha ~ lifespan$Freq);
summary(m);

#Plot regression analysis
plot(lifespan$RLS_Del_alpha ~ lifespan$Freq, ylab="Replicative Lifespan", xlab="Protein Interactions" );
abline( m, col="red");


###MORPHOLOGICAL PLASTICITY VS LIFESPAN REGRESSION ANALYSIS
#Match lifespan to mutant data
lifespan$stddev = scmd$stddev[match(lifespan$ORF, scmd$name)];

#Perform linear regression of lifespan to morphological plasticity and summarize results
m = lm(lifespan$RLS_Del_alpha ~ lifespan$stddev);
summary(m);

#Plot of regresssion analysis
plot(scmd$RLS_Del_alpha ~ scmd[,j]), ylab="Replicative Lifespan", xlab="Morphological Plasticity Robustness" );
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