###################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson & Nancy Yu
# Modified last: January 2020 by Annelie Johansson
###################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

library(biomaRt)
library(OCplus)
library(dplyr)

set.seed(97)
options(stringsAsFactors = FALSE)

#### Load and prepare data ----
load("input/20170404patient_array.RData")
patient_data <- select(ord_pat1, AgendiaRUNN, MammaPrint_Result, PAM50,
                       ERstatus, PRstatus, Ki67status, HER2status, grade, size20)
microarray <- ord_dat1
annotation <- annot1

dim(patient_data) # 652 x 9
allpat_info <- select(patient_data, AgendiaRUNN, ERstatus, MammaPrint_Result)

# All patients classidies as ultralow risk are ER-positive
all(allpat_info$ERstatus[which(allpat_info$MammaPrint_Result == "Ultra Low Risk")] == "Positive") # TRUE

# Select only ER-positive patients
erpat <- allpat_info[allpat_info$ERstatus == "Positive",]
erpat <- erpat[order(erpat$AgendiaRUNN),]
nrow(erpat) # 538
all(erpat$ERstatus == "Positive") # TRUE

table(erpat$MammaPrint_Result == "Ultra Low Risk")
# FALSE  TRUE 
# 440    98 

# Put in order: First utralow risk patients (UL), followed by all others (ER-positive, not UL)
# We want to compare UL vs ER-pos (not UL)
erpat_use <- erpat[c(which(erpat$MammaPrint_Result == "Ultra Low Risk"), which(erpat$MammaPrint_Result != "Ultra Low Risk")),]
all(erpat_use$MammaPrint_Result[1:98] == "Ultra Low Risk") # TRUE 

# Match gene expression data with patient data
dim(microarray) # 32146 x 652

all(erpat_use$AgendiaRUNN %in% colnames(microarray)) # TRUE
erdat <- microarray[, erpat_use$AgendiaRUNN]
identical(colnames(erdat), as.character(erpat_use$AgendiaRUNN)) # TRUE
dim(erdat) # 32146 x 538

# Check if all probes are unique? - yes!
nrow(erdat) # 32146
length(unique(rownames(erdat))) # 32146

### Filter on probes ----
# Remove probes with median expression less than the bottom 5% of the probes
# ie. Remove probes that are pretty much not expressed
med <- apply(erdat, 1, median, na.rm=T)
lqr <- quantile(erdat, na.rm=T, prob=0.05)
# plot(density(med))
# abline(v=lqr, col="red")
ind1 <- med > lqr
erdat_med <- erdat[ind1,]
nrow(erdat_med) # 30907
# 32146 - 30907 = 1239 probes removed

# Select probes with highest 75% IQR
# i.e. Remove probes with low variations
iqrange <- apply(erdat_med, 1, function(x) { IQR(x, na.rm=T) })
# plot(density(iqrange))
# abline(v=quantile(iqrange, prob=0.25) , col="red")
ind2 <- iqrange > quantile(iqrange, prob=0.25) 
erdat_iqr <- erdat_med[ind2, ]
nrow(erdat_iqr) # 23180
# 30907 - 23180 = 7727 probes removed

# Remove probes with NA values
erdat_complete <- erdat_iqr[complete.cases(erdat_iqr),]
nrow(erdat_complete) # 22394
# 23180 - 22394 = 786 probes removed 

# Match annotation data with gene expression data
all(rownames(erdat_complete) %in% rownames(annotation)) # TRUE
annot_complete <- annotation[rownames(erdat_complete), ]
nrow(annot_complete) # 22394

# Remove probes without gene names (from gene expression data + annotation data)
probes_without_genes <- rownames(annot_complete)[which(annot_complete$GeneName == "")]
length(probes_without_genes) # 1329
erdat_genes <- erdat_complete[-which(rownames(erdat_complete) %in% probes_without_genes),]
annot_genes <- annot_complete[-which(rownames(annot_complete) %in% probes_without_genes),]
nrow(erdat_genes) # 21065
nrow(annot_genes) # 21065
identical(rownames(erdat_genes), rownames(annot_genes)) # TRUE

## Remove AG_xxxx probes
probes_AG <- rownames(annot_genes)[grep("AG_", annot_genes$ProbeID)]
length(probes_AG) # 762
erdat_use <- erdat_genes[-which(rownames(erdat_genes) %in% probes_AG), ]
annot_use <- annot_genes[-which(rownames(annot_genes) %in% probes_AG), ]
nrow(erdat_use) # 20303
nrow(annot_use) # 20303
identical(rownames(erdat_use), rownames(annot_use)) # TRUE

### How many probes per gene? ----
probes_per_gene <- table(table(as.character(annot_use$GeneName)))
probes_per_gene
#     1     2     3     4     5     6     7     8     9    23 
# 11040  3250   652   127    36     7     3     3     1     1 
sum(probes_per_gene[-1]) # 4080 genes with multiple probes

### Run Differential Gene Analysis ----
# using function EOC in R package OCplus 
# OCplus is a package that uses t-statistics and False Discovery Rates (FDR) as cut offs
# for determining genes that are expressed at different levels between 2 groups

identical(colnames(erdat_use), as.character(erpat_use$AgendiaRUNN)) # TRUE
identical(rownames(erdat_use), rownames(annot_use)) # TRUE

# Perform differential analysis
erdat_eoc <- EOC(erdat_use, as.numeric(erpat_use$MammaPrint_Result == "Ultra Low Risk"))
erdat_eoc$GeneName <- annot_use$GeneName
erdat_eoc$ProbeID <- annot_use$ProbeID

# Find top genes, cutoff 0.001 FDR
erdat_eoc_topDE <- topDE(erdat_eoc, co=0.001)

# Using threshold of FDR < 0.001 results in n=793 significant probes (in n=706 genes)
nrow(erdat_eoc_topDE) # 793
length(unique(erdat_eoc_topDE$ProbeID)) # 793 probes
length(unique(erdat_eoc_topDE$GeneName)) # 706 genes

write.table(erdat_eoc_topDE, file="output/stable4.txt", sep="\t", quote=F, row.names=FALSE)
