library(RSparta)
library(limma)

sumlist <- read.maxq.summaryfile()
sumtab <- sumlist$sumtab
mq <- read.maxq.proteins(proteinGroups.filename = "proteinGroups.txt",log2.transform = T,null.zeros = T, is.label.free = T)

ptab <- mq$prot.table

ptab.nms <- mq$prot.colnms
ptab.metanms <- names(ptab)[! names(ptab) %in% ptab.nms]

ptab[ptab$Gene.names == 'Psen1',]
ptab.unnorm <- ptab
ptab.norm <-  normalizeBetweenArrays(as.matrix(ptab[,ptab.nms]),method="cyclicloess" )
ptab[,ptab.nms] <- ptab.norm[,ptab.nms]

ptab$AD2_3Mo_12fAD_Tube5 <- rowMeans(ptab[,c('AD2_3Mo_12fAD_Tube5_4p5ug','AD2_3Mo_12fAD_Tube5_2ug')],na.rm=T)
ptab.nms <- c('AD1_3Mo_14cF_Tube4_6ug', 'AD1_5Mo_11aC_Tube13',  'AD1_9Mo_16aD_Tube7',
              'AD2_3Mo_12fAD_Tube5',    'AD2_5Mo_11aABC_Tube14','AD2_9Mo_15bD_Tube8',
              'AD3_3Mo_14cA_Tube6',     'AD3_5Mo_12aC_Tube15',  'AD3_9Mo_16cCF_Tube9',
              'WT1_3Mo_14cC_Tube10',    'WT1_5Mo_11aD_Tube1',   'WT1_9Mo_16aA_Tube16',
              'WT2_3Mo_12fBE_Tube11',   'WT2_5Mo_11aA_Tube2',   'WT2_9Mo_16aF_Tube17',
              'WT3_3Mo_14cD_Tube12',    'WT3_5Mo_11aF_Tube3',   'WT3_9Mo_13BE_Tube18')
 
# remove one of the loading techreps
#ptab.nms <- ptab.nms[!ptab.nms %in%   c('AD2_3Mo_12fAD_Tube5_4p5ug')  ]                        
#ptab.nms <- ptab.nms[!ptab.nms %in%   c('AD2_3Mo_12fAD_Tube5_2ug')  ]                        

ptab$numvalid <- apply(ptab[,ptab.nms],1, function(x) sum(!is.na(x)) )
ptab$Amean <- apply(ptab[,ptab.nms],1, mean, na.rm=T)

ptab <- ptab[ptab$numvalid > 12 & ptab$Amean > 19.5,]
ptab[ptab$Gene.names == 'Psen1',]

gsub("(^[^_]+?)_.*","\\1",ptab.nms)

treatsib <-  gsub("(^[^_]+?)_.*","\\1",ptab.nms)
treat <-  factor(gsub("([^\\d])\\d+?","\\1",treatsib), levels=c('WT','AD'))



age <- factor(gsub(".*_(\\d+Mo).*","\\1",ptab.nms))

mouse <- gsub(".*_(\\d+.*)_Tube.*","\\1",ptab.nms)

geno <- factor(paste(treat,age,sep='.') , levels=c('WT.3Mo','AD.3Mo','WT.5Mo','AD.5Mo','WT.9Mo','AD.9Mo'))
#grps <- geno
#contrasts(geno) <- cbind(Time=c(-1,0,0,0,1,0),AD3mo=c(-1,1,0,0,0,0),AD5mo=c(0,0,-1,1,0,0),AD9mo=c(0,0,0,0,-1,1))
design <- model.matrix(~0 + geno)
# design <- model.matrix(~0+geno)
 colnames(design) <- levels(geno)
 
 rownames(design) <- ptab.nms
 
fit <- lmFit(ptab[,ptab.nms], design)
cont.wt <- makeContrasts("AD.9Mo-WT.9Mo", "AD.5Mo-WT.5Mo", "AD.3Mo-WT.3Mo", levels=design)
fit2 <- contrasts.fit(fit, cont.wt)

fit2 <- eBayes(fit2, trend = TRUE)

metanms <- c("uniprot.id","Gene.names")

topTableF(fit2, adjust="BH",genelist = ptab[,metanms],number=50)

ptopf <- topTableF(fit2, adjust="BH",genelist = ptab[,"Gene.names"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val


ftest.tab <- topTableF(fit2, adjust="BH",genelist = ptab[,metanms],number=Inf,p.value=1)

tab1 <- topTable(fit2, adjust="BH",genelist = ptab[,metanms], coef=1,p.value=1, number=Inf)
tab2 <- topTable(fit2, adjust="BH",genelist = ptab[,metanms], coef=2,p.value=1, number=Inf)
tab3 <- topTable(fit2, adjust="BH",genelist = ptab[,metanms], coef=3,p.value=1, number=Inf)
names(tab1)[3:8] <- paste(names(tab1)[3:8],"9Mo.ADvsWT",sep=".")
names(tab2)[3:8] <- paste(names(tab2)[3:8],"5Mo.ADvsWT",sep=".")
names(tab3)[3:8] <- paste(names(tab3)[3:8],"3Mo.ADvsWT",sep=".")

nms2merge <- c("uniprot.id","logFC","P.Value","adj.P.Val")

full.tab <- merge(ftest.tab,tab3[,c("uniprot.id","logFC.3Mo.ADvsWT","P.Value.3Mo.ADvsWT","adj.P.Val.3Mo.ADvsWT")], by = "uniprot.id",all = T, sort=F)
full.tab <- merge(full.tab,tab2[,c("uniprot.id","logFC.5Mo.ADvsWT","P.Value.5Mo.ADvsWT","adj.P.Val.5Mo.ADvsWT")], by = "uniprot.id",all = T, sort=F)
full.tab <- merge(full.tab,tab1[,c("uniprot.id","logFC.9Mo.ADvsWT","P.Value.9Mo.ADvsWT","adj.P.Val.9Mo.ADvsWT")], by = "uniprot.id",all = T, sort=F)
full.tab <- full.tab[order(full.tab$P.Value,decreasing=F),]

write.csv(full.tab,file="KB_MouseAD_empBayes_separateAgeandFtest.csv",row.names = F)

##


topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=c(2,5,6),number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=3,number=50)
topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef="AD.3Mo-WT.3Mo",number=50)

results <- decideTests(fit2)
vennDiagram(results)


design <- model.matrix(~treat*age)
rownames(design) <- ptab.nms
colnames(design) <- c('Intercept','treatAD','5Mo','9Mo','ADx5Mo','ADx9Mo')

fit <- lmFit(ptab[,ptab.nms], design)
#cont.wt <- makeContrasts("treatAD", levels=design)
#fit <- contrasts.fit(fit, cont.wt)

fit2 <- eBayes(fit,trend=T)
topTableF(fit2, adjust="BH",genelist = ptab[,"Gene.names"],number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=3,number=50)
topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef="AD.3Mo-WT.3Mo",number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=2,number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=2:6,p.value=0.2,number=Inf)
ptopf2 <- topTable(fit2, adjust="BH",genelist = ptab[,c("Gene.names","uniprot.acc")], coef=c(2,5,6),p.value=0.2,number=Inf)

##
design <- model.matrix(~age+treat)
rownames(design) <- ptab.nms
fit <- lmFit(ptab[,ptab.nms], design)
#cont.wt <- makeContrasts("treatAD", levels=design)
#fit <- contrasts.fit(fit, cont.wt)

fit2 <- eBayes(fit)
qqt(fit2$t,df=fit2$df.prior+fit2$df.residual,pch=16,cex=0.2)
abline(0,1)

topTableF(fit2, adjust.method = "BH",genelist = ptab[,"Gene.names"],number=50)
topTable(fit2, adjust.method = "BH",genelist = ptab[,"Gene.names"], coef = "treatAD",number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=3,number=50)
topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef="AD.3Mo-WT.3Mo",number=50)

topTable(fit2, adjust="BH",genelist = ptab[,"Gene.names"], coef=2,number=50)



#cont.matrix <- cbind(trt.ADvsWT=c(0,1,0,0,0,0),,Mu.SvsU=c(0,0,1,1),Diff=c(0,0,0,1))
#fit2 <- contrasts.fit(fit, cont.matrix)

#######
#######
treatments <- factor(c(1,1,2,2,3,3,4,4),labels=c("e10","E10","e48","E48"))
contrasts(treatments) <- cbind(Time=c(0,0,1,1),E10=c(0,1,0,0),E48=c(0,0,0,1))
design <- model.matrix(~treatments)


> colnames(design) <- c("Intercept","Time","E10","E48")