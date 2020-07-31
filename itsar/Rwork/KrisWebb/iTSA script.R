library("DEP")
library("dplyr")
library("EnrichmentBrowser")
library("limma")


data <- read.table(file = "proteinGroups_filt.txt",
                   header = TRUE,sep = "\t",quote = "\"'",dec = ".",numerals = c("warn.loss"),
                   row.names = NULL,na.strings = c("NA","NaN","Infinite"))                           #Read proteinGroups file

#data <- filter(data, Reverse != "+", Potential.contaminant != "+") #filter for contaminant proteins and decoy database hits, which are indicated by
#data <- filter(data, Reverse != "+")                                                                   #"+" in the columns "Potential.contaminants" and "Reverse", respectively.

data$Gene.names %>% duplicated() %>% any()                         # Are there any duplicated gene names?
#any(duplicated(data$Gene.names))


data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>%   # Make a table of duplicated gene names
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Generate a SummarizedExperiment object using an experimental design
TMT_columns <- grep("Reporter.intensity", colnames(data_unique)) # get TMT column numbers

experimental_design <-read.csv(file = "ExperimentalDesign.csv", header = TRUE)

experimental_design$label <- as.character(experimental_design$label)             #correct matrix data type for make_se
experimental_design$condition <- as.character(experimental_design$condition)
experimental_design$replicate <- as.numeric(experimental_design$replicate)

str(experimental_design)

data_se <- make_se(data_unique, TMT_columns, experimental_design)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)
#data_norm <- data_filt
#assay(data_norm) <- normalizeBetweenArrays(assay(data_filt), method = "cyclicloess", cyclic.method="affy")
#assay(data_norm) <- normalizeBetweenArrays(assay(data_filt), method = "quantile")

plot_normalization(data_filt, data_norm)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp <- impute(data_norm, fun = "knn", rowmax = 0.9)

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "C")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0.2))

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0.98, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 3, col_limit = 1, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 1.5, show_row_names = FALSE)

# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "D_vs_C", label_size = 3, add_names = TRUE)

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

# Generate a wide data.frame
df_wide <- get_df_wide(dep)

# Generate a long data.frame
df_long <- get_df_long(dep)

# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")

# These data can be loaded in future R sessions using this command
load("data.RData")

#Write resuts to CSV
write.csv(df_wide, file="df_wide.csv")