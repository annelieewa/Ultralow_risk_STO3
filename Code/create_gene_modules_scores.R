#################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson & Nancy Yu -- modified code from Nick Tobin
# Modified last: January 2020 by Annelie Johansson
#################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

library("biomaRt")
library("genefu")
options(stringsAsFactors = FALSE)

#### Load data -----
load("input/20170404patient_array.RData")
microarray <- ord_dat1
annotation <- annot1

# Load all_modules: a list of 24 gene modules with the probe IDs, Entrez IDs, and coefficients
load("input/all_modules.RData")

#### Add Entrez ID:s to annot1 using R package biomaRt ----
# we need one column with "EntrezGene.ID" to calculate gene signature scores

# Select biomart database
hg19 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

# Retreive annotations from Ensembl, this takes some time
hg19_allhgnc <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","entrezgene", "description", "gene_biotype", "refseq_mrna"), 
                      filters = "with_hgnc", values = TRUE, mart = hg19)
nrow(hg19_allhgnc) # 76158
# remove NA:s for entrezgene
hg19_entrez_NAs_removed  <- hg19_allhgnc[!is.na(hg19_allhgnc$entrezgene), ]
nrow(hg19_entrez_NAs_removed) # 64061

# Obtain newer Entrez IDs from Ensembl.
# Not all microarray probes have annotated Entrez IDs.
# Some IDs might be outdated.
annotation$EntrezGene.ID <- hg19_entrez_NAs_removed[match(annotation$GeneName, hg19_entrez_NAs_removed$hgnc_symbol), ]$entrezgene
# Some Refseq IDs in annot1 are obsolete
length(which(is.na(annotation$EntrezGene.ID))) # 4504 - i.e. 14%
length(which(!is.na(annotation$EntrezGene.ID))) # 27642

#### Calculate gene modules scores using genefu ----
# https://www.bioconductor.org/packages/release/bioc/html/genefu.html
# https://rdrr.io/bioc/genefu/
# https://www.ncbi.nlm.nih.gov/pubmed/26607490

## Function that uses genefu package's sig.score to calculate score for a given gene module
mod_score2 <- function(eSet, annotation, moduleOb){
  require(genefu)
  require(Biobase)
  mod_ob = sig.score(
    x = moduleOb,
    data = t(eSet),
    annot = annotation,
    do.mapping = TRUE,
    signed = TRUE,
    verbose = TRUE)        
  return(data.frame(mod_ob$score))
}

## Function to convert module scores to tertile numbers 0, 1, or 2 
tertiles_to_integer = function (dfOb){  
  ifelse((dfOb < 0) == TRUE, 
         (ifelse(dfOb < quantile(dfOb, probs = 1/3), 0, ifelse(dfOb >= quantile(dfOb, probs = 1/3) & dfOb < quantile(dfOb, probs = 2/3), 1, 2))),
         (ifelse(dfOb < quantile(dfOb, probs = 1/3), 0, ifelse(dfOb >= quantile(dfOb, probs = 1/3) & dfOb < quantile(dfOb, probs = 2/3), 1, 2))))
}

## Initialize genefu_scores data frame and calculate 
genefu_scores <- data.frame(matrix(NA, nrow = ncol(microarray), ncol = 24), row.names = colnames(microarray))
for(i in 1:24) {
  genefu_scores[i] <- mod_score2(eSet = microarray, annotation = annotation, moduleOb = all_modules[[i]])
}

## Number of Entrez IDs within each gene module gene list that Genefu recognizes
# Genefu output:
# probe candidates: 51/56   91%   Gene70
# probe candidates: 70/72   97%   CIN70
# probe candidates: 50/50   100%  STROMA1
# probe candidates: 64/68   94%   STROMA2
# probe candidates: 6/6     100%  IMMUNE1
# probe candidates: 91/95   96%   IMMUNE2
# probe candidates: 208/276 75%   RAS 
# probe candidates: 359/379 95%   MAPK
# probe candidates: 174/205 85%   PTEN
# probe candidates: 399/503 79%   AKT_MTOR
# probe candidates: 683/712 96%   IGF1
# probe candidates: 45/47   96%   SRC
# probe candidates: 157/191 82%   MYC
# probe candidates: 192/222 86%   E2F3
# probe candidates: 55/65   85%   BETAC
# probe candidates: 102/128 80%   GGI
# probe candidates: 225/278 81%   PIK3CA
# probe candidates: 451/469 96%   ESR1
# probe candidates: 27/28   96%   ERBB2
# probe candidates: 224/229 98%   AURKA
# probe candidates: 64/68   94%   PLAU
# probe candidates: 14/14   100%  VEGF
# probe candidates: 91/95   96%   STAT1
# probe candidates: 10/10   100%  CASP3

colnames(genefu_scores) <- names(all_modules)
module_variance <- apply(genefu_scores, 2, sd) 
summary(module_variance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07353 0.12818 0.18029 0.18857 0.24683 0.50346

#### Divide gene module scores into tertiles ----
# Categorizes the three tertiles to two groups – the worst vs the other two
# (which depends upon if high or low score is associated with bad or poor prognosis) 

tertiles_all <- data.frame(matrix(NA, nrow =  nrow(genefu_scores), ncol =  ncol(genefu_scores)), row.names = colnames(microarray))

# Loop this function over each column  
for (i in 1: ncol(genefu_scores)){  
  tertiles_all[i] <- tertiles_to_integer(genefu_scores[,i])
}
names(tertiles_all) = names(genefu_scores)
num_tertile <- apply(tertiles_all, 2, function(x) table(x))
# check so that we have equal numbers in each group
num_tertile

# Order the modules alphabetically
genefu_scores <- genefu_scores[ ,order(colnames(genefu_scores))]
tertiles_all <- tertiles_all[ ,order(colnames(tertiles_all))]  
identical(colnames(genefu_scores), colnames(tertiles_all)) # TRUE

colnames(tertiles_all)

## We want to compare the worst tertile vs the other 2

# Gene modules were high scores are associated with poor prognosis (= Highest tertile is worst)
High_genes <- c("AKT_MTOR", "AURKA", "BETAC", "CIN70", "E2F3", "ERBB2", "GENE70", "GGI", "IGF1", "IMMUNE1",
                "IMMUNE2", "MAPK", "MYC", "PLAU", "PTEN", "RAS", "SRC", "STROMA1", "STROMA2")
# Gene modules were low scores are associated with poor prognosis (= Lowest tertile is worst)
Low_genes <- c("CASP3", "ESR1", "PIK3CA", "STAT1", "VEGF")

# Convert tertiles to binary number (worst tertile vs. the other 2 combined)
tertiles_hi <- tertiles_all[, High_genes]
tertiles_low <- tertiles_all[, Low_genes]
tertiles_hi_bin <- apply(tertiles_hi, 2, function(x) { ifelse(x == 2, 2, 1) })
tertiles_low_bin <- apply(tertiles_low, 2, function(x) { ifelse(x == 0, 1, 2) }) 
tertiles_bin <- cbind(tertiles_hi_bin, tertiles_low_bin)
tertiles_bin <- tertiles_bin[,order(colnames(tertiles_bin))]

dim(genefu_scores) # 652 x 24
dim(tertiles_all) # 652 x 24
dim(tertiles_bin) # 652 x 24

head(genefu_scores)
head(tertiles_all)
head(tertiles_bin) 

## tertiles_bin is the data we are using in for downstream analysis
## remove gene modules not analysed
# PLAU = same as STROMA2
# STAT1 = same as IMMUNE2
# not including: CIN70, GENE70, GGI
gene_modules <- tertiles_bin[, -which(colnames(tertiles_bin) %in% c("PLAU", "STAT1", "CIN70",
                                                                    "GENE70", "GGI"))]

save(gene_modules, file = "output/gene_modules.RData")
