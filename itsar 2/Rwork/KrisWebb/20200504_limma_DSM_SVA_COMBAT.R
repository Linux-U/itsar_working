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

proteinGroups <- dplyr::filter(proteinGroups, Reverse != "+", Potential.contaminant != "+", Only.identified.by.site != "+")



#Filter for protines with 2 or more unique peptides

proteinGroups <- dplyr::filter(proteinGroups, Unique.peptides.O > 1)
proteinGroups <- dplyr::filter(proteinGroups, Unique.peptides.R > 1)


#Filter for protines with 2 or more reporter ion reads

proteinGroups <- dplyr::filter(proteinGroups, Reporter.intensity.count.1.O > 2)
proteinGroups <- dplyr::filter(proteinGroups, Reporter.intensity.count.1.R > 2)


## Make row ID unique name

# Are there any duplicated gene names?

proteinGroups$Gene.names %>% duplicated() %>% any() 



# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.

proteinGroups <- make_unique(proteinGroups, "Gene.names", "Protein.IDs", delim = ";")

rownames(proteinGroups) <- proteinGroups[,"name"]

# divide columns by control 1 to make between plex comparisons

proteinGroups <- transform(proteinGroups, Reporter.intensity.2.O = Reporter.intensity.2.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.3.O = Reporter.intensity.3.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.4.O = Reporter.intensity.4.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.5.O = Reporter.intensity.5.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.6.O = Reporter.intensity.6.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.7.O = Reporter.intensity.7.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.8.O = Reporter.intensity.8.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.9.O = Reporter.intensity.9.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.10.O = Reporter.intensity.10.O / Reporter.intensity.1.O)

proteinGroups <- transform(proteinGroups, Reporter.intensity.11.O = Reporter.intensity.11.O / Reporter.intensity.1.O)



proteinGroups <- transform(proteinGroups, Reporter.intensity.2.R = Reporter.intensity.2.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.3.R = Reporter.intensity.3.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.4.R = Reporter.intensity.4.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.5.R = Reporter.intensity.5.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.6.R = Reporter.intensity.6.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.7.R = Reporter.intensity.7.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.8.R = Reporter.intensity.8.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.9.R = Reporter.intensity.9.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.10.R = Reporter.intensity.10.R / Reporter.intensity.1.R)

proteinGroups <- transform(proteinGroups, Reporter.intensity.11.R = Reporter.intensity.11.R / Reporter.intensity.1.R)










#Build an ExpressionSet from proteinGroups data.frame



#1 Assay data



#create a minimal ExpressionSet object with just the assay data wanted

#assaydata <- data.matrix(proteinGroups[,82:103])

assaydata <- data.matrix(proteinGroups[, c(83:92  , 94:103)])

DSM2020_R <- ExpressionSet(assayData=assaydata)



#2 Phenotypic data.  Data describing treatment, control, batch, and other covariates.



#Make a phenotype table describing your experiment.

phenotable <- matrix(c("DH","C","DH","C","DH","C","DH","C","DH","C","DL","C","DL","C","DL","C","DL","C","DL","C","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2"),nrow = 20, ncol = 2)

colnames(phenotable) <- c("Treatment","Batch")

rownames(phenotable) <- colnames(exprs(DSM2020_R))

phenotable <- as.data.frame(phenotable)

phenotable <- new("AnnotatedDataFrame", data = phenotable)



# Verify row names of phenotable are the same as the column names of the assay data of the expression set.

all(rownames(phenotable)==colnames(exprs(DSM2020_R)))

DSM2020_R <- ExpressionSet(assayData=assaydata, phenoData=phenotable)



#3 Feature data.

featuretable <- proteinGroups[c(1:81,104:147)]

featuretable <- as.data.frame(featuretable)

featuretable <- new("AnnotatedDataFrame", data = featuretable)



# Verify row names of featuretable are the same as the row names of the assay data of the expression set.

all(rownames(featuretable)==rownames(exprs(DSM2020_R)))

DSM2020_R <- ExpressionSet(assayData=assaydata, phenoData=phenotable, featureData=featuretable)

#ExpressionSet is ready



#Load dataset

pdata = pData(DSM2020_R)
edata = exprs(DSM2020_R)
fdata = fData(DSM2020_R)



#Log2 transform all edata

edata = log2(edata)



#Quantile normalization

norm_edata = normalize.quantiles(as.matrix(edata))
colnames(norm_edata) <- colnames(edata)
row.names(norm_edata) <- row.names(edata)



# By default calculates the euclidean distance between rows.

dist1 = dist(t(edata))



#Heatmap of euclidean distance between samples.

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)



#hierarchical clustering of samples.

hclust1 = hclust(dist1)
plot(hclust1,hang=-1)



#Adjusting for batch effects with Combat

batch = pdata$Batch
modcombat = model.matrix(~1, data=pdata)
modcancer = model.matrix(~Treatment, data=pdata)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

dist1 = dist(t(combat_edata))
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)
hclust1 = hclust(dist1)
plot(hclust1,hang=-1)



#Move back to protein groups due to errors in SVA.  Hopefuly this hack it to work

proteinGroups[,83:92] <- combat_edata[,1:10]
proteinGroups[,94:103] <- combat_edata[,11:20]




#Build an ExpressionSet from proteinGroups data.frame



#1 Assay data



#create a minimal ExpressionSet object with just the assay data wanted

assaydata2 <- data.matrix(proteinGroups[,c(83,85,87,89,91,94,96,98,100,102)])

DSM2020_O_vs_R <- ExpressionSet(assayData=assaydata2)



#2 Phenotypic data.  Data describing treatment, control, batch, and other covariates.



#Make a phenotype table describing your experiment.

phenotable2 <- matrix(c("DH","DH","DH","DH","DH","DL","DL","DL","DL","DL"),nrow = 10, ncol = 1)

colnames(phenotable2) <- c("Treatment")

rownames(phenotable2) <- colnames(exprs(DSM2020_O_vs_R))

phenotable2 <- as.data.frame(phenotable2)

phenotable2 <- new("AnnotatedDataFrame", data = phenotable2)



# Verify row names of phenotable are the same as the column names of the assay data of the expression set.

all(rownames(phenotable2)==colnames(exprs(DSM2020_O_vs_R)))

DSM2020_O_vs_R <- ExpressionSet(assayData=assaydata2, phenoData=phenotable2)



#3 Feature data.

featuretable2 <- proteinGroups[c(1:81,104:147)]

featuretable2 <- as.data.frame(featuretable2)

featuretable2 <- new("AnnotatedDataFrame", data = featuretable2)



# Verify row names of featuretable are the same as the row names of the assay data of the expression set.

all(rownames(featuretable2)==rownames(exprs(DSM2020_O_vs_R)))

DSM2020_O_vs_R <- ExpressionSet(assayData=assaydata2, phenoData=phenotable2, featureData=featuretable2)

#ExpressionSet is ready



#Load dataset

pdata2 = pData(DSM2020_O_vs_R)
edata2 = exprs(DSM2020_O_vs_R)
fdata2 = fData(DSM2020_O_vs_R)


# By default calculates the euclidean distance between rows.

dist1 = dist(t(edata2))



#Heatmap of euclidean distance between samples.

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)



#hierarchical clustering of samples.

hclust1 = hclust(dist1)
plot(hclust1,hang=-1)



#Remove unknown batch variables with surrogate variable analysis.

mod = model.matrix(~Treatment,data=pdata2)
mod0 = model.matrix(~1, data=pdata2)
sva1 = sva(edata2,mod,mod0)



#Add the surrogate variables to the model matrix and perform the model fit.

modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata2))



#statistics with limma

fit_limma = lmFit (edata2, modsv)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma)

names(ebayes_limma)
hist(ebayes_limma$p.value[,"TreatmentDL"],col=2)


 
#adjust P-values for multiple hypothesis testing.  Benjamini Hochberg

qstats_obj <- p.adjust(ebayes_limma$p.value[,"TreatmentDL"], method = "BH", length(ebayes_limma$p.value[,"TreatmentDL"]))
sum(qstats_obj < 0.05)




#Stop work 20200525


#Make a proteingroups object with limma p values, q values, and expression diffrences

expdiff <- as.matrix(rowMeans(edata2[,c(1,2,3,4,5)])-rowMeans(edata2[,c(6,7,8,9,10)]))
colnames(expdiff) <- "log2 fold change"

p.value <- ebayes_limma$p.value[,"TreatmentDL"]
p.value <- as.data.frame(p.value)
names(p.value)[1] <- "p value vehicle"

qstats_obj <- as.matrix(qstats_obj)
qstats_obj <- as.data.frame(qstats_obj)
names(qstats_obj)[1] <- "q value"

negLogq = -log2(qstats_obj)
negLogq <- as.data.frame(negLogq)
names(negLogq)[1] <- "-log2(q)"


stats_table <- cbind( expdiff, p.value, qstats_obj, negLogq)

#proteinGroups[,93:103] <- norm_edata[,1:11]
proteinGroups_limma <- merge(proteinGroups, stats_table, by=0)

write.csv(proteinGroups_limma, file ="proteinGroups_limma_SVA_O_vs_R.csv")







