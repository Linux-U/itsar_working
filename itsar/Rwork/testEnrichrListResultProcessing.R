
detach("package:RSparta", unload = T)

library(RSparta)
#data("enrichr.libraries")

library(enrichR)
# selected.databases <- c("Allen_Brain_Atlas_down","Allen_Brain_Atlas_up",
#                         "GO_Biological_Process_2017","GO_Molecular_Function_2017", "GO_Cellular_Component_2017", 
#                         "Reactome_2016","WikiPathways_2016","KEGG_2016","HMDB_Metabolites","Humancyc_2016","huMAP",
#                         "Jensen_COMPARTMENTS","Jensen_DISEASES","Jensen_TISSUES",
#                         "BioPlex_2017","Genome_Browser_PWMs","OMIM_Disease",
#                         "CORUM","LINCS_L1000_Chem_Pert_down"  ,                     
#                         "LINCS_L1000_Chem_Pert_up",  "LINCS_L1000_Kinase_Perturbations_down"  ,          
#                         "LINCS_L1000_Kinase_Perturbations_up" , "LINCS_L1000_Ligand_Perturbations_down" ,"LINCS_L1000_Ligand_Perturbations_up",
#                         "MSigDB_Oncogenic_Signatures","MSigDB_Computational",
#                         "ENCODE_Histone_Modifications_2015","ENCODE_TF_ChIP-seq_2015","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
#                         "Epigenomics_Roadmap_HM_ChIP-seq","Chromosome_Location","ChEA_2016",
#                         "Genes_Associated_with_NIH_Grants",
#                         "Disease_Perturbations_from_GEO_down","Disease_Perturbations_from_GEO_up",
#                         "Disease_Signatures_from_GEO_down_2014","Disease_Signatures_from_GEO_up_2014",
#                         "Drug_Perturbations_from_GEO_down","Drug_Perturbations_from_GEO_up",
#                         "DrugMatrix",
#                         "GeneSigDB","NCI-60_Cancer_Cell_Lines","NCI-Nature_2016","OMIM_Expanded",
#                         "GTEx_Tissue_Sample_Gene_Expression_Profiles_down","GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
#                         "NURSA_Human_Endogenous_Complexome",
#                         "Kinase_Perturbations_from_GEO_up","Kinase_Perturbations_from_GEO_down",
#                         "Ligand_Perturbations_from_GEO_up","Ligand_Perturbations_from_GEO_down",
#                         "PPI_Hub_Proteins","SILAC_Phosphoproteomics",
#                         "Panther_2016","Pfam_InterPro_Domains","Phosphatase_Substrates_from_DEPOD",
#                         "Single_Gene_Perturbations_from_GEO_down","Single_Gene_Perturbations_from_GEO_up","TargetScan_microRNA",
#                         "TF-LOF_Expression_from_GEO","Transcription_Factor_PPIs",
#                         "BioCarta_2016",
#                         "Cancer_Cell_Line_Encyclopedia","RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
#                         "TRANSFAC_and_JASPAR_PWMs")
# 
# selected.databases <- unique(selected.databases)
# elibs <- listEnrichrDbs()
# sort(elibs$libraryName)
# 
# selected.databases[!selected.databases %in% elibs$libraryName]
# 
# selected.databases[!selected.databases %in% names(enrichr.libraries)]
# 
# 
# sort(elibs$libraryName[!elibs$libraryName %in% selected.databases ])
# 
# libs2add <- download.enrichr.libraries(libraryNames =  selected.databases[!selected.databases %in% names(enrichr.libraries)] )
libs2add <- download.enrichr.libraries(libraryNames = c("KEA_2015","ESCAPE"))

# enrichr.libraries <- c(enrichr.libraries,libs2add)
# enrichr.libraries <- enrichr.libraries[unique(sort(names(enrichr.libraries)))]
# save(enrichr.libraries, file='D:/Will/R/ENRICHR_DATABASES/enrichr_libraries_download_best.RData')
# enrichr.libraries <- libraries
#devtools::use_data(enrichr.libraries, pkg = "D:/Will/R/Functions/RSparta", internal = F, compress = "bzip2")

#load('D:/Will/R/ENRICHR_DATABASES/enrichr_libraries_download_best.RData')

print.enriched.set <- function(enrich.list, colinds = c(1,2,3,4,5,6), maxchar = 35) {
  for(x in 1:length(names(enrich.list))) {
    readline(prompt = names(enrich.list)[x])
    
    if(sum(enrich.list[[x]]$adj.P.Val < 0.05) > 0) {
      mydf <- enrich.list[[x]]
      mydf$Genes <- substr(as.character(mydf$Genes),1,maxchar)
      maxpr <- min(sum(mydf$adj.P.Val < 0.05), 10)
      print(mydf[1:maxpr,colinds])  
    }
    
  }
  
}

orggenes <- tnew$Gene.name[!is.na(tnew$Gene.name) ]

hupgenes <- tnew$Gene.name[!is.na(tnew$Gene.name) & !is.na(tnew$adj.P.Val.harm) & !is.na(tnew$logFC.harm) & tnew$adj.P.Val.harm < 0.1 & tnew$logFC.harm > 0]

hdown <- tnew$Gene.name[!is.na(tnew$Gene.name) & !is.na(tnew$adj.P.Val.harm) & !is.na(tnew$logFC.harm) & tnew$adj.P.Val.harm < 0.1 & tnew$logFC.harm < 0]

T21highgenes <- tnew$Gene.name[!is.na(tnew$Gene.name) & !is.na(tnew$adj.P.Val.cellh50) & !is.na(tnew$logFC.cellh50) & tnew$adj.P.Val.cellh50 < 0.1 & tnew$logFC.cellh50 < 0]

T21lowgenes <- tnew$Gene.name[!is.na(tnew$Gene.name) & !is.na(tnew$adj.P.Val.cellh50) & !is.na(tnew$logFC.cellh50) & tnew$adj.P.Val.cellh50 < 0.1 & tnew$logFC.cellh50 > 0]
#try the intersection:
harmup.norm <- intersect(T21lowgenes, hupgenes)
harmdn.norm <- intersect(T21highgenes, hdown)

length(hupgenes)
length(hdown)
length(T21highgenes)
length(T21lowgenes)
writeLines(hupgenes, con = 'hupgenes_20170810.txt')
writeLines(hdown, con = 'hdown_20170810.txt')
writeLines(T21highgenes, con = 'T21highgenes_20170810.txt')
writeLines(T21lowgenes, con = 'T21lowgenes_20170810.txt')

testenrich <- RSparta::perform.fisher.enrich( termlist = enrichr.libraries$GO_Biological_Process_2017, 
                                              diffgeneset = T21highgenes, 
                                              background = orggenes, minOverlapSize = 4, minTermListSize = 8, maxTermListSize = 800)



t21htest  <-  plyr::llply(enrichr.libraries["RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"], .progress = "text",
                          perform.fisher.enrich, 
                          diffgeneset = T21highgenes, 
                          background = orggenes, minOverlapSize = 8, minTermListSize = 10, maxTermListSize = 800)

minOverlapSize <- 8
minTermListSize <- 10
maxTermListSize <- 1800

rm(T21.CO.Hup0.1.enrichr.fulldflist,T21.CO.Hdown0.1.enrichr.fulldflist,T21.CO.T21high.0.1.enrichr.fulldflist,T21.CO.T21low0.1.enrichr.fulldflist,
   T21.CO.HupNORM0.1.enrichr.fulldflist,T21.CO.HdownNORM0.1.enrichr.fulldflist)


T21.CO.Hup0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, hupgenes, orggenes, minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize) #,  .progress = "text")

T21.CO.Hdown0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, hdown, orggenes,  minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize)

T21.CO.T21high.0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, T21highgenes, orggenes,  minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize)

T21.CO.T21low0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, T21lowgenes, orggenes, minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize)

T21.CO.HupNORM0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, harmup.norm, orggenes,  minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize)

T21.CO.HdownNORM0.1.enrichr.fulldflist <-  plyr::llply(enrichr.libraries, perform.fisher.enrich, harmdn.norm, orggenes,  minOverlapSize = minOverlapSize, minTermListSize = minTermListSize, maxTermListSize = maxTermListSize)

print.enriched.set(T21.CO.Hup0.1.enrichr.fulldflist)

print.enriched.set(T21.CO.T21low0.1.enrichr.fulldflist)

print.enriched.set(T21.CO.T21high.0.1.enrichr.fulldflist)

print.enriched.set(T21.CO.Hdown0.1.enrichr.fulldflist)




print.enriched.set(T21.CO.HupNORM0.1.enrichr.fulldflist)


print.enriched.set(T21.CO.HdownNORM0.1.enrichr.fulldflist)





###
library(dplyr)

# # Because the character columns in the data frames in these lists
# # were incorrectly created as factors, need to convert. 
# # TO DO: FIX the perform.fisher.enrich() function to create them as character columns
#  
# T21.CO.Hup0.1.enrichr.fulldflist <- lapply( T21.CO.Hup0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 
# T21.CO.Hdown0.1.enrichr.fulldflist <- lapply( T21.CO.Hdown0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 
# T21.CO.T21high.0.1.enrichr.fulldflist <- lapply( T21.CO.T21high.0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 
# T21.CO.T21low0.1.enrichr.fulldflist <- lapply( T21.CO.T21low0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 
# T21.CO.HupNORM0.1.enrichr.fulldflist <- lapply( T21.CO.HupNORM0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 
# T21.CO.HdownNORM0.1.enrichr.fulldflist <- lapply( T21.CO.HdownNORM0.1.enrichr.fulldflist, function(xdf)    xdf %>% mutate_if(is.factor, as.character)  )
# 

save(hupgenes,hdown,T21highgenes,T21lowgenes,
     prots, tnew, tnew.nonimpute,
     T21.CO.Hup0.1.enrichr.fulldflist,
     T21.CO.T21low0.1.enrichr.fulldflist,
     T21.CO.T21high.0.1.enrichr.fulldflist,
     T21.CO.Hdown0.1.enrichr.fulldflist, T21.CO.HupNORM0.1.enrichr.fulldflist, 
     T21.CO.HdownNORM0.1.enrichr.fulldflist,file = 'TMB_Paper_EnrichRanalysis_20170812.RData'  )

print.enriched.set(T21.CO.Hup0.1.enrichr.fulldflist, maxchar = 35)

T21.CO.Hup.table <- dplyr::bind_rows(T21.CO.Hup0.1.enrichr.fulldflist, .id = "enrichrlib")

# Note on categories/libraries to keep
# RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO T21high, harmine down T21 high: 1137        Trisomy_21_Primary_Fibroblasts_GSE55504_up 1.169176e-10 1.484853e-07      51      496 TAGLN;FBN1;CSPG4;ACTC1;CPZ;FLNC;DCN
# 
# Trisomy_21_Primary_Fibroblasts_GSE55504_up
# WikiPathways_2016
# Disease_Perturbations_from_GEO_down
# Disease_Perturbations_from_GEO_up
# Drug_Perturbations_from_GEO_2014
# Drug_Perturbations_from_GEO_down
# Drug_Perturbations_from_GEO_up
# ENCODE_TF_ChIP-seq_2015
# ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
# Jensen_COMPARTMENTS
# Jensen_DISEASES **Down_syndrome 3.096190e-03 0.0154809484       8       89 COL6A1;GART;CRELD1;PFKL;SOD1;HMGN1;
# Jensen_TISSUES
# KEGG_2016
# LINCS_L1000_Ligand_Perturbations_up
#Ligand_Perturbations_from_GEO_up
#Ligand_Perturbations_from_GEO_down
# Kinase_Perturbations_from_GEO_down (BRD4 inhibition, FGFR1 activation)
#Kinase_Perturbations_from_GEO_up (LRRK2 active mutant, TGFBR2 KO, CDK8 KD, MAP3K7 KD)
#LINCS_L1000_Chem_Pert_down
# NCI-60_Cancer_Cell_Lines 
# MSigDB_Computational
# Panther_2016
#Reactome_2016
# Humancyc_2016
# RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO
# Single_Gene_Perturbations_from_GEO_down
# Single_Gene_Perturbations_from_GEO_up
# Chromosome_Location
# LINCS_L1000_Chem_Pert_down 1756     LJP006_LNCAP_24H-CHIR-99021-10 2.947862e-07 0.0002925604      13       95 TAGLN;DCN;MGP;RARRES2;COL6A1;MDK;S1
# LINCS_L1000_Chem_Pert_up
# LINCS_L1000**
# # BioPlex_2017
# huMAP 
# GO_*s
# Maybes
# Ligand_Perturbations_from_GEO_down
# Ligand_Perturbations_from_GEO_up
# NCI-Nature_2016
# Epigenomics_Roadmap_HM_ChIP-seq
# 
# 

# bigdf <- dplyr::filter(bigdf, !enrichrlib %in% c("GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
#                                                  "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",
#                                                  "Jensen_COMPARTMENTS",
#                                                  "ENCODE_Histone_Modifications_2015",
#                                                  "ENCODE_Histone_Modifications_2015",
#                                                  "Epigenomics_Roadmap_HM_ChIP-seq",
#                                                  "GeneSigDB","huMAP","PPI_Hub_Proteins","Genes_Associated_with_NIH_Grants",
#                                                  "Cancer_Cell_Line_Encyclopedia"))

conv.list.to.df <- function(listofdfs) {
  bigdf <- dplyr::bind_rows(listofdfs, .id = "enrichrlib")
  bigdf$global.adj.P.Val <- p.adjust(bigdf$adj.P.Val)
  bigdf <- bigdf %>% dplyr::arrange(P.Value)
  return(bigdf)
}

T21.CO.Hup.table <- conv.list.to.df(T21.CO.Hup0.1.enrichr.fulldflist)

T21.CO.Hdn.table <- conv.list.to.df(T21.CO.Hdown0.1.enrichr.fulldflist)

T21.CO.T21.HIGH.table <- conv.list.to.df(T21.CO.T21high.0.1.enrichr.fulldflist)

T21.CO.T21.LOW.table <- conv.list.to.df(T21.CO.T21low0.1.enrichr.fulldflist)

T21.CO.Hup.T21lowNORM.table <- conv.list.to.df(T21.CO.HupNORM0.1.enrichr.fulldflist)

T21.CO.Hdn.T21highNORM.table <- conv.list.to.df(T21.CO.HdownNORM0.1.enrichr.fulldflist)

dflist <- list(Harmine.UP = T21.CO.Hup.table,
               Harmine.DN = T21.CO.Hdn.table,
               T21.HIGH = T21.CO.T21.HIGH.table,
               T21.LOW = T21.CO.T21.LOW.table)


topcats <- sort(unique(unlist(lapply(dflist, function(xdf) unique(xdf$enrichrlib[xdf$global.adj.P.Val < 0.05])))))

topdf <- dplyr::bind_rows(dflist, .id = "DEX.comp.name") %>% dplyr::filter(global.adj.P.Val < 0.1)

topdf <- topdf %>% arrange(enrichrlib, Term, DEX.comp.name)


#                
# > sort(names(enrichr.libraries))
# [1] "Achilles_fitness_decrease"                        "Achilles_fitness_increase"                       
# [3] "Allen_Brain_Atlas_down"                           "Allen_Brain_Atlas_up"                            
# [5] "BioPlex_2017"                                     "ChEA_2016"                                       
# [7] "Chromosome_Location"                              "CORUM"                                           
# [9] "Disease_Perturbations_from_GEO_down"              "Disease_Perturbations_from_GEO_up"               
# [11] "Disease_Signatures_from_GEO_down_2014"            "Disease_Signatures_from_GEO_up_2014"             
# [13] "Drug_Perturbations_from_GEO_2014"                 "Drug_Perturbations_from_GEO_up"                  
# [15] "DrugMatrix"                                       "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"       
# [17] "ENCODE_Histone_Modifications_2015"                "ENCODE_TF_ChIP-seq_2015"                         
# [19] "Epigenomics_Roadmap_HM_ChIP-seq"                  "Genes_Associated_with_NIH_Grants"                
# [21] "GeneSigDB"                                        "GO_Biological_Process_2017"                      
# [23] "GO_Cellular_Component_2017"                       "GO_Molecular_Function_2017"                      
# [25] "GTEx_Tissue_Sample_Gene_Expression_Profiles_down" "GTEx_Tissue_Sample_Gene_Expression_Profiles_up"  
# [27] "HMDB_Metabolites"                                 "Humancyc_2016"                                   
# [29] "huMAP"                                            "Jensen_COMPARTMENTS"                             
# [31] "Jensen_DISEASES"                                  "Jensen_TISSUES"                                  
# [33] "KEGG_2016"                                        "Kinase_Perturbations_from_GEO_down"              
# [35] "Kinase_Perturbations_from_GEO_up"                 "Ligand_Perturbations_from_GEO_down"              
# [37] "Ligand_Perturbations_from_GEO_up"                 "MSigDB_Computational"                            
# [39] "MSigDB_Oncogenic_Signatures"                      "NCI-60_Cancer_Cell_Lines"                        
# [41] "NCI-Nature_2016"                                  "NURSA_Human_Endogenous_Complexome"               
# [43] "OMIM_Expanded"                                    "Panther_2016"                                    
# [45] "Pfam_InterPro_Domains"                            "Phosphatase_Substrates_from_DEPOD"               
# [47] "PPI_Hub_Proteins"                                 "Reactome_2016"                                   
# [49] "SILAC_Phosphoproteomics"                          "Single_Gene_Perturbations_from_GEO_down"         
# [51] "Single_Gene_Perturbations_from_GEO_up"            "TargetScan_microRNA"                             
# [53] "TF-LOF_Expression_from_GEO"                       "Transcription_Factor_PPIs"                       
# [55] "Virus_Perturbations_from_GEO_up"                  "Virus_Perturbations_from_GEO_up"                 
# [57] "WikiPathways_2016"    














