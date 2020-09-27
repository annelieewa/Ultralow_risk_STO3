##### Project: Ultralow risk STO3 ###############################################
# fishers_test.R
# Author: Annelie Johansson & Adina Iftimi
# Modified last: January 2020 by Annelie Johansson
#################################################################################

library(dplyr)

#### Load data ----

# Load data file containing: annot1, ord_dat1, ord_pat1
load("20170404patient_array.RData")

# Load gene module scores
load("Genefu_modulescores.RData")

# Select only clinical patient data that is needed
ord_pat2 <- dplyr::select(ord_pat1, AgendiaRUNN, MammaPrint_Result, PAM50, ERstatus, PRstatus, HER2status, Ki67status, size20, grade)
dim(ord_pat2) # 652 x 9

# Select only ER-positive patients
ord_pat2 <- subset(ord_pat2, ERstatus=="Positive")
nrow(ord_pat2) # 538

## Reference group: Ultralow risk patients (UL)
ord_pat2$UL <- ord_pat2$MammaPrint_Result == "Ultra Low Risk"

## ER-positive group: ER-positive, not UL 
ord_pat2$ERpos <- (ord_pat2$MammaPrint_Result != "Ultra Low Risk")

## Luminal A group: Luminal A subtype, not UL
ord_pat2$LumA <- (ord_pat2$PAM50 == "LumA" & ord_pat2$MammaPrint_Result != "Ultra Low Risk")

## Luminal B group: Luminal B subtype, not UL
ord_pat2$LumB <- (ord_pat2$PAM50 == "LumB" & ord_pat2$MammaPrint_Result != "Ultra Low Risk")

pat_use <- ord_pat2
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

# Change some variables names and categories (order from best to worse prognosis-marker)
pat_use$PR <- pat_use$PRstatus 
pat_use$HER2 <- pat_use$HER2status 
pat_use$Ki67 <- pat_use$Ki67status 
pat_use$Size <- pat_use$size20 
pat_use$Grade <- pat_use$grade 

pat_use$PR <- factor(pat_use$PR, levels=c("Positive", "Negative"), labels=c("positive", "negative"))
pat_use$HER2 <- factor(pat_use$HER2, levels=c("Negative", "Positive"), labels=c("negative", "positive"))
pat_use$Ki67 <- factor(pat_use$Ki67, levels=c("Negative", "Positive"), labels=c("low", "high"))
pat_use$Size <- factor(pat_use$Size, levels=c("0", "1"), labels=c("<=20 mm", ">20 mm"))
pat_use$Grade <- factor(pat_use$Grade, levels=c("1", "2", "3"), labels=c("1", "2", "3"))

#### Add gene module scores ----

# Match with patients data
tertiles_bin <- tertiles_bin[match(pat_use$AgendiaRUNN, rownames(tertiles_bin)),]
identical(pat_use$AgendiaRUNN, rownames(tertiles_bin)) # TRUE

# Change names on categories
tertiles_bin[which(tertiles_bin == "1")] <- "low"
tertiles_bin[which(tertiles_bin == "2")] <- "high"

# Get correct variable names and categories (order from best to worse prognosis-marker)
tertiles_bin <- as.data.frame(tertiles_bin)

# Gene modules were low scores are associated with poor prognosis (= Lowest tertile is worst)
Low_genes <- c("CASP3", "ESR1", "IMMUNE2", "PIK3CA", "STAT1", "VEGF")

for(i in 1:ncol(tertiles_bin)){
  if(colnames(tertiles_bin)[i] %in% Low_genes){
    tertiles_bin[,i] <- factor(tertiles_bin[,i], levels = c("high", "low"))
  }else{ tertiles_bin[,i] <- factor(tertiles_bin[,i], levels = c("low", "high")) }
}


# Add to pat_use
tertiles_bin$AgendiaRUNN <- rownames(tertiles_bin)
pat_use <- merge(pat_use, tertiles_bin, by = "AgendiaRUNN", all.x = TRUE)

pat_use <- subset(pat_use, )
variables <- c("PR", "HER2", "Ki67", "Size", "Grade", colnames(genefu_scores))


####  Create table for fisher's test and mosaic plot ----
# We want to compare:
#  - Ultralow patients (UL) vs ER-positive patients (not UL)
#  - Ultralow patients vs Luminal A patients (ER-positive, not UL)
#  - Ultralow patients vs Luminal B patients (ER-positive, not UL)
table_mosaic <- NULL

for(v in 1:length(variables)){
  var_use <- variables[v]
  data <- pat_use
  #Rvar <- data[, grep(var_use, colnames(data))]
  var <- data[, which(colnames(data) == var_use)]
  
  if(any(is.na(var))){ data <- data[-which(is.na(var)), ] }
  var <- data[, which(colnames(data) == var_use)]
  
  #Rvar <- as.factor(data[, grep(var_use, colnames(data))])
  nr_table <- NULL
  for(i in 1:length(levels(var))){
    nr <- c(length(which(var == levels(var)[i] & data$UL == TRUE)),
            length(which(var == levels(var)[i] & data$ERpos == TRUE)),
            length(which(var == levels(var)[i] & data$LumA == TRUE)),
            length(which(var == levels(var)[i] & data$LumB == TRUE)) )
    nr_table <- rbind(nr_table, nr)
  }
  temp <- cbind.data.frame(cbind.data.frame(rep(var_use, length(levels(var))), levels(var)), nr_table)
  colnames(temp) <- c("Characteristics", "level", "UL", "ERp", "LumA", "LumB")
  table_mosaic <- rbind(table_mosaic, temp)
}

#### Create table with results from Fisher's test ----

list_pvals <- vector("list", length(variables))

set.seed(123)
for(v in 1:length(variables)){
  tempERp <- table_mosaic[which(table_mosaic$Characteristics == variables[v]), c("UL", "ERp")]
  tempLumA <- table_mosaic[which(table_mosaic$Characteristics == variables[v]), c("UL", "LumA")]
  tempLumB <- table_mosaic[which(table_mosaic$Characteristics == variables[v]), c("UL", "LumB")]
  pvals <- c(fisher.test(tempERp, alternative = "two.sided")$p.value,
             fisher.test(tempLumA, alternative = "two.sided")$p.value,
             fisher.test(tempLumB, alternative = "two.sided")$p.value)
  list_pvals[[v]] <- pvals
  names(list_pvals)[v] <- variables[v]
}

table_pvals <- matrix(unlist(list_pvals), ncol = 3, byrow=TRUE)
rownames(table_pvals) <- names(list_pvals)
colnames(table_pvals) <- c("UL_vs_ERp", "UL_vs_LumA", "UL_vs_LumB")

# write.table(table_mosaic, "table_mosaic.txt", row.names = FALSE)
# write.table(table_pvals, "table_fishers_test.txt", row.names = TRUE)

