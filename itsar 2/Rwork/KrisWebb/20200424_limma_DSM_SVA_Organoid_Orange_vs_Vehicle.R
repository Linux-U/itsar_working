library("DEP")
library("dplyr")
library("EnrichmentBrowser")
library("limma")
library("Biobase")
library("dendextend")
library("RSkittleBrewer")
library("preprocessCore")
library("vsn")
library("broom")
library("snpStats")
library("sva")
library("genefilter")




#Read proteinGroups file

proteinGroups <- read.table(file = "proteinGroups.txt",
                   header = TRUE,sep = "\t",quote = "\"'",dec = ".",numerals = c("warn.loss"),
                   row.names = NULL,na.strings = c("NA","NaN","Infinite"))   

proteinGroups <- as.data.frame(proteinGroups)

#Filter dataset

#Filter for contaminant proteins, reverse hits, and only identifed my site.

proteinGroups <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")



#Filter for protines with 2 or more unique peptides

proteinGroups <- filter(proteinGroups, Unique.peptides.O > 1)



#Filter for protines with 2 or more unique peptides

proteinGroups <- filter(proteinGroups, Reporter.intensity.count.1.O > 2)



## Make row ID unique name

# Are there any duplicated gene names?

proteinGroups$Gene.names %>% duplicated() %>% any() 



# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.

proteinGroups <- make_unique(proteinGroups, "Gene.names", "Protein.IDs", delim = ";")

rownames(proteinGroups) <- proteinGroups[,"name"]



#Build an ExpressionSet from proteinGroups data.frame



#1 Assay data



#create a minimal ExpressionSet object with just the assay data wanted

assaydata <- data.matrix(proteinGroups[,82:92])

DSM2020_O <- ExpressionSet(assayData=assaydata)



#2 Phenotypic data.  Data describing treatment, control, batch, and other covariates.



#Make a phenotype table describing your experiment.

phenotable <- matrix(c("C","D","C","D","C","D","C","D","C","D","C"))

colnames(phenotable) <- c("Treatment")

rownames(phenotable) <- colnames(exprs(DSM2020_O))

phenotable <- as.data.frame(phenotable)

phenotable <- new("AnnotatedDataFrame", data = phenotable)



# Verify row names of phenotable are the same as the column names of the assay data of the expression set.

all(rownames(phenotable)==colnames(exprs(DSM2020_O)))

DSM2020_O <- ExpressionSet(assayData=assaydata, phenoData=phenotable)



#3 Feature data.

featuretable <- proteinGroups[c(1:26,104:114,126,128:147)]

featuretable <- as.data.frame(featuretable)

featuretable <- new("AnnotatedDataFrame", data = featuretable)



# Verify row names of featuretable are the same as the row names of the assay data of the expression set.

all(rownames(featuretable)==rownames(exprs(DSM2020_O)))

DSM2020_O <- ExpressionSet(assayData=assaydata, phenoData=phenotable, featureData=featuretable)

#ExpressionSet is ready



#Exploratory Analysis

trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)



#Load dataset
pdata=pData(DSM2020_O)
edata=exprs(DSM2020_O)
fdata = fData(DSM2020_O)



#boxplot to look at log2 distribution of samples

boxplot(log2(edata+1),col=2)



#histograms

par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)

par(mfrow=c(1,1))



#density plot

plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)
lines(density(log2(edata[,3]+1)),col=4)
lines(density(log2(edata[,4]+1)),col=5)
lines(density(log2(edata[,5]+1)),col=6)
lines(density(log2(edata[,6]+1)),col=7)
lines(density(log2(edata[,7]+1)),col=8)
lines(density(log2(edata[,8]+1)),col=9)
lines(density(log2(edata[,9]+1)),col=10)
lines(density(log2(edata[,10]+1)),col=11)
lines(density(log2(edata[,11]+1)),col=12)



#qq-plot. Compare distributions of measurements before normalization.

qqplot(log2(edata[,1]+1), log2(edata[,8]+1),col=3)
abline(c(0,1))



#Bland Altman plot

mm = log2(edata[,1]+1) - log2(edata[,4]+1)
aa = log2(edata[,1]+1) + log2(edata[,4]+1)
plot(aa,mm,col=2)



#Heatmap

ematrix = as.matrix(edata)
heatmap(ematrix)

#change up the colors

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp)

#No clustering

heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)



#Clustering



#Log2 transform all edata

edata = log2(edata)



# By default calculates the euclidean distance between rows.

dist1 = dist(t(edata))



#Heatmap of euclidean distance between samples.

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)



#hierarchical clustering of samples.

hclust1 = hclust(dist1)
plot(hclust1,hang=-1)

#color the dendrogram

dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,4,col=1:4)
plot(dend)



#Kmeans clustering

kmeans1 = kmeans(edata,centers=3)
names(kmeans1)

#look at the cluster centers

matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)

#how many belong to each cluster

table(kmeans1$cluster)

#cluster the data together and plot

heatmap(as.matrix(edata)[order(kmeans1$cluster),],col=colramp,Colv=NA,Rowv=NA)



#Dimension reduction

#Mean center rowdata

edata_centered = edata - rowMeans(edata)



#calculate singular vectors

svd1 = svd(edata_centered)

names(svd1)



#Plot d singular values

plot(svd1$d,ylab="Singular value",col=2)



#Plot Percent Variance Explained

plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)



#Plot top two singular vector decomposition values

par(mfrow=c(1,2))
plot(svd1$v[,1],col=2,ylab="1st PC")
plot(svd1$v[,2],col=2,ylab="2nd PC")



#Plot SV1 vs. SV2

par(mfrow=c(1,1))

plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC",xlab="1st PC")

plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$Treatment))


#Boxplot of SV1 and Treatment

boxplot(svd1$v[,1] ~ pdata$Treatment,border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$Treatment)),col=as.numeric(pdata$Treatment))



#principal component vs singular vector.  

pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])

#If you subtract the column means insted of the row means they are the same

edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2)



#Normalization 

#Plot density distribution of samples

plot(density(edata[,1]),col=colramp[7],lwd=5,ylim=c(0,.2))
for(i in 2:11){lines(density(edata[,i]),lwd=5,col=colramp[i])}



#Quantile normalization

norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=5,ylim=c(0,.20))
for(i in 2:11){lines(density(norm_edata[,i]),lwd=5,col=colramp[i+2])}
colnames(norm_edata) <- colnames(edata)
row.names(norm_edata) <- row.names(edata)

#VSN normalization

#vsnnorm_edata = justvsn(as.matrix(2^edata))
#plot(density(vsnnorm_edata[,1]),col=colramp[1],lwd=5,ylim=c(0,0.2))
#for(i in 2:11){lines(density(vsnnorm_edata[,i]),lwd=5,col=colramp[i+2])}


#Matching distributions leaves variability. Note normalization will not remove batch effects.

#plot(norm_edata[1,],col=as.numeric(pdata$Treatment))
#plot(vsnnorm_edata[1,],col=as.numeric(pdata$Treatment))

svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata$Treatment), pch = 19, cex=2)



#Remove unknown batch variables with surrogate variable analysis.

mod = model.matrix(~Treatment,data=pdata)
mod0 = model.matrix(~1, data=pdata)
sva1 = sva(norm_edata,mod,mod0)




#Add the surrogate variables to the model matrix and perform the model fit.

modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(norm_edata))


#T-statistic

tstats_obj = rowttests(norm_edata,pdata$Treatment)
names(tstats_obj)
hist(tstats_obj$statistic,col=2)
hist(tstats_obj$p.value,col=2)

#statistics with limma

fit_limma = lmFit (norm_edata, modsv)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma)

names(ebayes_limma)
hist(ebayes_limma$statistic,col=2)
hist(ebayes_limma$p.value[,"TreatmentD"],col=2)




plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)



#adjust P-values

qstats_obj <- p.adjust(ebayes_limma$p.value[,"TreatmentD"], method = "BH", length(ebayes_limma$p.value[,"TreatmentD"]))
sum(qstats_obj < 0.05)



#Adjusted p-values from limma.

#limma_pvals_adj = topTable(ebayes_limma,number=dim(edata)[1])$adj.P.Val
#hist(limma_pvals_adj,col=2)
#quantile(limma_pvals_adj)

#Make a proteingroups object with limma p values, q values, and expression diffrences

expdiff <- as.matrix(rowMeans(norm_edata[,c(1,3,5,7,9,11)])-rowMeans(norm_edata[,c(2,4,6,8,10)]))
colnames(expdiff) <- "log2 fold change"

p.value <- ebayes_limma$p.value[,"TreatmentD"]
names(p.value)[1] <- "p value"

qstats_obj <- as.matrix(qstats_obj)
names(qstats_obj)[1] <- "q stat"

negLogq = -log2(qstats_obj)
names(negLogq)[1] <- "-log2 q"


stats_table <- cbind( expdiff, p.value, qstats_obj, negLogq)

proteinGroups[,82:92] <- norm_edata[,1:11]
proteinGroups_limma <- merge(proteinGroups, stats_table, by=0)

write.csv(proteinGroups_limma, file ="proteinGroups_limma_vehicle_vs_1uM_bCX.csv")
