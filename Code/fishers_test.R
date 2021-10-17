###################################################################################
# Manuscript: Clinical and Molecular Characteristics of ER-Positive Ultralow Risk
#            Breast Cancer Tumors Identified by the 70-Gene Signature
# Author: Annelie Johansson / Adina Iftimi
# Modified last: August 2021 by Annelie Johansson
###################################################################################

setwd("/Volumes/Annelie Encrypted/Projects/Ultralow_risk_genes/Analysis/")

library(dplyr)

## Load data ----
load("input/20170404patient_array.RData")
patient_data <- select(ord_pat1, AgendiaRUNN, MammaPrint_Result, PAM50,
                       ERstatus, PRstatus, Ki67status, HER2status, grade, size20)
microarray <- ord_dat1
annotation <- annot1

## Load gene modules 
load("output/gene_modules.RData")

# All patients classifies as ultralow risk are ER-positive
all(patient_data$ERstatus[which(patient_data$MammaPrint_Result == "Ultra Low Risk")] == "Positive") # TRUE

## Select only ER-positive patients
identical(colnames(microarray), patient_data$AgendiaRUNN) # TRUE
patient_data <- subset(patient_data, ERstatus=="Positive")
nrow(patient_data) # 538
microarray <- microarray[,patient_data$AgendiaRUNN]
identical(colnames(microarray), patient_data$AgendiaRUNN) # TRUE

## Reference group: Only Ultralow patients ####
patient_data$UL <- patient_data$MammaPrint_Result == "Ultra Low Risk"

## ER-positive group: ER-positive, not UL ####
patient_data$ERpos <- (patient_data$MammaPrint_Result != "Ultra Low Risk")

## Luminal A group: Luminal A, not UL ####
patient_data$LumA <- (patient_data$PAM50 == "LumA" & patient_data$MammaPrint_Result != "Ultra Low Risk")

## Luminal B group: Luminal B, not UL ####
patient_data$LumB <- (patient_data$PAM50 == "LumB" & patient_data$MammaPrint_Result != "Ultra Low Risk")

pat_use <- patient_data
nrow(pat_use) # 538
table(pat_use$UL, pat_use$ERpos)
#           FALSE TRUE
# FALSE     0  440
# TRUE     98    0
table(pat_use$UL, pat_use$LumA)
#         FALSE TRUE
# FALSE   191  249
# TRUE     98    0
table(pat_use$UL, pat_use$LumB)
#         FALSE TRUE
# FALSE   318  122
# TRUE     98    0

## Fix variables
pat_use$size20 <- as.character(pat_use$size20)
pat_use$size20[which(is.na(pat_use$size20))] <- "Unknown"
pat_use$size20 <- factor(pat_use$size20, levels = c("0", "1", "Unknown"))
table(pat_use$size20, useNA = "always")
#   0       1 Unknown    <NA> 
# 437      95       6       0 

pat_use$grade <- as.character(pat_use$grade)
pat_use$grade[which(is.na(pat_use$grade))] <- "Unknown"
pat_use$grade <- factor(pat_use$grade, levels = c("1", "2", "3", "Unknown"))
table(pat_use$grade, useNA = "always")
#    1       2       3 Unknown    <NA> 
#  116     338      76       8       0 

pat_use$PRstatus <- factor(pat_use$PRstatus, levels = c("Positive", "Negative", "Unknown"))
table(pat_use$PRstatus, useNA = "always")
# Positive Negative  Unknown     <NA> 
#    367      162        9        0 

pat_use$HER2status <- factor(pat_use$HER2status, levels = c("Negative", "Positive", "Unknown"))
table(pat_use$HER2status, useNA = "always")
# Negative Positive  Unknown     <NA> 
#    513       24        1        0

pat_use$Ki67status <- factor(pat_use$Ki67status, levels = c("Negative", "Positive", "Unknown"))
table(pat_use$Ki67status, useNA = "always")
# Negative Positive  Unknown     <NA> 
#     395      117       26        0 

pat_use$PAM50 <- factor(pat_use$PAM50, levels = c("LumA", "LumB", "Normal", "Her2", "Basal"))
table(pat_use$PAM50, useNA = "always")
# LumA   LumB Normal   Her2  Basal   <NA> 
# 336    126     48     21      7      0 


## Add gene modules

# match with patients data
gene_modules <- gene_modules[match(pat_use$AgendiaRUNN, rownames(gene_modules)),]
identical(pat_use$AgendiaRUNN, rownames(gene_modules)) # TRUE
gene_modules <- as.data.frame(gene_modules)

# add to pat_use
gene_modules$AgendiaRUNN <- rownames(gene_modules)
pat_use <- merge(pat_use, gene_modules, by="AgendiaRUNN", all.x = TRUE)

#### FISHERS TEST ####

gene_modules_names <- colnames(gene_modules)[-20]
variables <- c("size20", "grade", "PRstatus", "HER2status", "Ki67status", "PAM50",
               gene_modules_names)

# Compare UL - ERpos (not UL)
# Compare UL - LumA (not UL)
# Compare UL - LumB (not UL)

table_fisher <- NULL

for(v in 1:length(variables)){
  var_use <- variables[v]
  data <- pat_use
  var <- data[, grep(var_use, colnames(data))]
  
  if(var_use %in% gene_modules_names){
    l1 <- length(which(var == 1))
    l2 <- length(which(var == 2))
    if(l1 > l2){
      var <- factor(var, levels = c(1, 2), labels = c("Low", "High"))
    }else{ var <- factor(var, levels = c(2, 1), labels = c("High", "Low")) }
  }
  
  # table with actual number 
  nr_table <- NULL
  for(i in 1:length(levels(var))){
    nr <- c(length(which(var == levels(var)[i] & data$UL == TRUE)),
            length(which(var == levels(var)[i] & data$ERpos == TRUE)),
            length(which(var == levels(var)[i] & data$LumA == TRUE)),
            length(which(var == levels(var)[i] & data$LumB == TRUE)) )
    nr_table <- rbind(nr_table, nr)
  }
  
  if(any(levels(var) == "Unknown")){
    temp <- nr_table[1:length(levels(var))-1,]
  }else{ temp <- nr_table }
  
  # percentages
  perc_table <- formatC(t(t(temp) / colSums(temp))*100, format = "f", digits = 1)
  if("Unknown" %in% levels(var)){ perc_table <- rbind(perc_table, rep("-", 4)) } 
  
  # fisher's test
  pvals <- c(fisher.test(temp[,c(1,2)], alternative = "two.sided")$p.value,
             fisher.test(temp[,c(1,3)], alternative = "two.sided")$p.value,
             fisher.test(temp[,c(1,4)], alternative = "two.sided")$p.value)
    pvals <- formatC(pvals, format = "f", digits = 3)
    pvals[which(pvals == "0.000")] <- "<0.001"
  
  # put together number + percentages + p-values
  tab <- cbind(paste0(nr_table[,1], "_(", perc_table[,1], ")"),
                paste0(nr_table[,2], "_(", perc_table[,2], ")"),
                c(pvals[1], rep("", length(levels(var))-1)),
                paste0(nr_table[,3], "_(", perc_table[,3], ")"),
                c(pvals[2], rep("", length(levels(var))-1)),
                paste0(nr_table[,4], "_(", perc_table[,4], ")"),
                c(pvals[3], rep("", length(levels(var))-1)))
  
  colnames(tab) <- c("UL", "ERpos", "P1", "LumA", "P2", "LumB", "P3")
  rownames(tab) <- levels(var)
  table_fisher <- rbind(table_fisher, c(var_use, rep("", 6)), tab)
}

nrow(table_fisher) # 84

table1_save <- table_fisher[1: which(table_fisher[,1] == gene_modules_names[1])-1,]
stable2_save <- table_fisher[which(table_fisher[,1] == gene_modules_names[1]):nrow(table_fisher),]

# ## Save tables
write.table(table1_save, "output/table1.txt", sep = "\t") #, row.names = FALSE)
write.table(stable2_save, "output/stable3.txt", sep = "\t") #, row.names = TRUE)


