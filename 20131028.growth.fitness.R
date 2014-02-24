rm(list=ls())
setwd("~/projects/0.network.aging.prj/0.ppi.reliability.simulation")
list.files(path="data/", pattern="Regression")
my.files = c("Regression_Tc1_het.txt", "Regression_Tc1_hom.txt", "Regression_Tc2_het.txt", "Regression_Tc2_hom.txt")
# I removed the single quotes in these files

tb1 = read.table(paste('data/',my.files[1],sep=''), header=T, sep='\t', fill=T)
tb2 = read.table(paste('data/',my.files[2],sep=''), header=T, sep='\t', fill=T)
tb3 = read.table(paste('data/',my.files[3],sep=''), header=T, sep='\t', fill=T)
tb4 = read.table(paste('data/',my.files[4],sep=''), header=T, sep='\t', fill=T)

tb1$AnnotationFlag = ifelse( is.na(tb1$YPD), 'questionable1', 'verified1' )
table(tb1$AnnotationFlag)
# questional verified 
# 174         5744 

tb3$AnnotationFlag = ifelse( is.na(tb3$YPD), 'questionable2', 'verified2' )
table(tb3$AnnotationFlag)
# quetionable    verified 
# 152        5766 
hist(tb3$YPD)

tb2$GrowthFlag = ifelse( is.na(tb2$YPD), 'nogrowth1', 'growth1' )
table(tb2$GrowthFlag)
# growth nogrowth 
# 4659     1259 

tb4$GrowthFlag = ifelse( is.na(tb4$YPD), 'nogrowth2', 'growth2' )
table(tb4$GrowthFlag)
# growth nogrowth 
# 4718     1200 

tb12 = merge(tb1,tb2, by.x='orf', by.y='orf')
tb12 = tb12[, c(1,2,grep('Flag',names(tb12)))]

tb34 = merge(tb3,tb4, by.x='orf', by.y='orf')
tb34 = tb34[, c(1,2,grep('Flag',names(tb34)))]

names(tb12) = c("orf",'gene','Anno1','Growth1')
names(tb34) = c('orf','gene','Anno2','Growth2')
str(tb12)
str(tb34)

tb = merge(tb12,tb34, by.x='orf', by.y='orf')
head(tb)
tb$name.check= ifelse(tb$gene.x == tb$gene.y, T, F)
table(tb$name.check)
# tb[ tb$name.check==F, ]
length( grep('question', tb$Anno1) ) #174
length( grep('question', tb$Anno2) ) #152

#tbReal = tb[ tb$Anno1=='verified1' & tb$Anno2=='verified2', ]
tbReal = tb[ tb$Anno1=='verified1' | tb$Anno2=='verified2', ]  #At least one Het growth. The other one might be an error?!
length( grep('question', tbReal$Anno1) ) #28
length( grep('question', tbReal$Anno2) ) #6

table(tbReal$Growth1, tbReal$Growth2)
#           growth2 nogrowth2
# growth1      4552         5
# nogrowth1      63      1152

# Verify some manually on SGD.
tbReal[tbReal$Growth1=='nogrowth1' & tbReal$Growth2=='nogrowth2',  ] [1:10,]
tbReal[tbReal$Growth1=='growth1' & tbReal$Growth2=='nogrowth2',  ] 

tbReal$essenflag = ifelse( tbReal$Growth1=='growth1' & tbReal$Growth2=='growth2', 'nonessential', 
                       ifelse(tbReal$Growth1=='nogrowth1' & tbReal$Growth2=='nogrowth2', 'essential', 'abnormal') )
head(tbReal)
table(tbReal$essenflag) 
# abnormal    essential nonessential 
# 68         1152         4552 
#Good, consistent results. Use these for further analysis. 2013 Oct 29

write.csv(tbReal, "SummaryRegressionHetHom2013Oct29.csv",  quote=F, row.names=F )

