mqdatacolnms <- c("null","rep1","rep2","rep3LF")


inpath <- "D:/Will/PROJECTS/INDUSTRY_COLLABORATIONS/SYROS_COLLABORATION/combined_allPhospho_HL60s_10-25-2017/HL60_Proteome/"
fname <- "proteinGroups_ProteomeOnly_AllReps_1-22-2018_PerseusAnnotated.txt"
fpath <- paste(inpath,fname,sep="")

protcols <- unlist(strsplit( readLines(fpath ,n  =  1),"\t"))

syrprots <- read.table(fpath,header=T,stringsAsFactors=F, na.strings="NaN",
                       quote = "",
                       sep="\t", fill = T)

syrprots.test <- read.maxq.proteins(proteinGroups.filename = fpath,
                                    datacolnames = mqdatacolnms,
                                    uniprot.acc.column = "Majority protein IDs",
                                    uniprot.id.column = "UniProt names",
                                    datacolumnprefix = "",
                                    log2.transform = F,
                                    is.label.free = F)

