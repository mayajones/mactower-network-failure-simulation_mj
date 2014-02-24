# use long andn not-long categories

library(e1071)
library(pcurve);

gasch = read.csv("gasch00.tab", sep="\t", header=T );
row.names(gasch) = gasch[,1]

sub.g = gasch[, seq(4,176)]
sub.g[ is.na(sub.g) ] = 0; #this means log(mut/wt)=0, i.e. no change of expression

row.names(sub.g) = gasch[,1]
gpca = pca(sub.g );
row.names(gpca$pcs) = gasch[,1]

#write.table(gpca$pcs, "gasch.pca.tab", col.name=F, sep="\t", quote=F);
gpcs = data.frame( gpca$pcs) ;
gpcs[ 1:5, 1:2] 
summary( lm( x[,1] ~ x[,2] ) ) #check the pca results

tb = read.csv("RLS.of.564.gene.deletion.BY.Managbanag08Plos.csv")
tb$ORF = as.character(tb$ORF)
row.names(tb) = as.character( tb$ORF )

tb$ratio = tb$RLS_Del_alpha / tb$BY4742 #R2 = 0.83

#cutoffs = quantile(tb$RLS_Del_alpha, prob=c(0.75, 0.90), na.rm=T);
#cutoffs = quantile(tb$ratio, prob=c(0.75, 0.90), na.rm=T); #not working
cutoffs = quantile(tb$ratio, prob=c(0.25, 0.90), na.rm=T); #not working

tb$score = NA;
tb$score[tb$ratio< cutoffs[1] ] = "notlong"
tb$score[tb$ratio> cutoffs[2] ] = "long"
table(tb$score)

shared.orfs = intersect(tb$ORF, row.names(gpcs)) #shared orfs bw RLS and Gasch datasets
sub = tb[shared.orfs, ]
#sub = tb[ ! is.na(tb$score), ]

###svm prediction
y = factor( sub$score )
#x = gpcs[ shared.orfs, 1:100]
x = gpcs[ shared.orfs, ]
dim(x)

model = svm( x, y);
     # test with train data
     pred <- predict(model, x)
     # (same as:) #     pred <- fitted(model)
     
     # Check accuracy:
     table( pred, y)  


#how many pca should we choose for prediction?
summary( lm( sub$ratio ~ x[,1] + x[,2] + x[,3] + x[,4] + x[,5] +x[,6]  ) )

pred.gasch = predict(model, gpcs);
table( pred.gasch) #




####################loose the different cutoff values

#decide to use 10% as cutoff 
tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.1 ] = "increase";
tb$score[ (tb$arls/26.5)<0.9 ] = "decrease";
factor(tb$score);
y = factor( tb$score )
model = svm( x, y);
pred <- predict(model, x)
table( pred, y)  #Ha, this is 100% correct.

 pred.gasch = predict(model, gpcs[,1:43]);

summary(pred.gasch); #shit 4573 decrease genes, only 5 increase genes from training set.


#### try again
tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.075 ] = "increase";  ##this give 26 increase from whole genome
tb$score[ (tb$arls/26.5)<0.85 ] = "decrease";
factor(tb$score);
y = factor( tb$score )
model = svm( x, y);
pred <- predict(model, x)
table( pred, y)  #Ha, this is 100% correct.

 pred.gasch = predict(model, gpcs[,1:43]);

summary(pred.gasch); 
n = row.names(gpcs);
plus.orf = n[ pred.gasch=="increase"]

gasch[plus.orf, 2]

write.table( tb, "kaeberlein.2.tab", quote=F, col.names=F, row.names=F,sep="\t");

write.table( gasch[plus.orf, 2], "111505.predicted.rls.extension.orfs.tab" , quote=F, row.name=F,col.names=F);



