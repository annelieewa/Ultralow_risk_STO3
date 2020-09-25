##### Project: Ultralow risk STO3 ###############################################
# eTable2_add_hallmark_GO_groups.R
# Author: Annelie Johansson & Christina Yau
# Modified last: January 2020 by Annelie Johansson
#################################################################################

library(sigPathway)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("sigPathway")

## Hallmark genes and GO biological processes genes from here:
# http://software.broadinstitute.org/gsea/msigdb/index.jsp

#### Load data ----

# Load list of genes associated with ultralow risk
diff_genes <- read.table("topgenes_seed97_793significant.txt", sep='\t', header=TRUE, stringsAsFactor=F)

# Read in Hallmarks gene sets
hallmarks <- gmtToG("Data/h.all.v6.2.symbols.gmt.txt") 

# Read in GO BP gene groups to be used
load("Data/GO_use.RData") # loads c5names_use and c5probes_use

nrow(diff_genes) # 793
length(unique(diff_genes$ProbeID)) # 793


#### START HALLMARKS ####
# Find which genes belongs to which Hallmarks

hnames <- sapply(hallmarks, function(x) {x$src})
hprobes <- sapply(hallmarks, function(x) {x$probes})
length(hprobes) # 50
length(hnames) # 50 

## Merge some of the hallmarks ----

# We put together ESTROGEN_RESPONSE_EARLY and ESTROGEN_RESPONSE_LATE
hnames[which(hnames == "HALLMARK_ESTROGEN_RESPONSE_EARLY")] <- "HALLMARK_ESTROGEN_RESPONSE"
hprobes[[which(hnames == "HALLMARK_ESTROGEN_RESPONSE")]] <- c(hprobes[[which(hnames == "HALLMARK_ESTROGEN_RESPONSE")]],
                                                              hprobes[[which(hnames == "HALLMARK_ESTROGEN_RESPONSE_LATE")]])
hprobes <- hprobes[-which(hnames == "HALLMARK_ESTROGEN_RESPONSE_LATE")]
hnames <- hnames[-which(hnames == "HALLMARK_ESTROGEN_RESPONSE_LATE")]

length(hprobes) # 49
length(hnames) # 49

# We put together KRAS_SIGNALING_DN and KRAS_SIGNALING_UP
hnames[which(hnames == "HALLMARK_KRAS_SIGNALING_DN")] <- "HALLMARK_KRAS_SIGNALING"
hprobes[[which(hnames == "HALLMARK_KRAS_SIGNALING")]] <- c(hprobes[[which(hnames == "HALLMARK_KRAS_SIGNALING")]],
                                                              hprobes[[which(hnames == "HALLMARK_KRAS_SIGNALING_UP")]])
hprobes <- hprobes[-which(hnames == "HALLMARK_KRAS_SIGNALING_UP")]
hnames <- hnames[-which(hnames == "HALLMARK_KRAS_SIGNALING_UP")]

length(hprobes) # 48
length(hnames) # 48

# We put together MYC_TARGETS_V1 and MYC_TARGETS_V2
hnames[which(hnames == "HALLMARK_MYC_TARGETS_V1")] <- "HALLMARK_MYC_TARGETS"
hprobes[[which(hnames == "HALLMARK_MYC_TARGETS")]] <- c(hprobes[[which(hnames == "HALLMARK_MYC_TARGETS")]],
                                                           hprobes[[which(hnames == "HALLMARK_MYC_TARGETS_V2")]])
hprobes <- hprobes[-which(hnames == "HALLMARK_MYC_TARGETS_V2")]
hnames <- hnames[-which(hnames == "HALLMARK_MYC_TARGETS_V2")]

length(hprobes) # 47
length(hnames) # 47

# We put together UV_RESPONSE_DN and UV_RESPONSE_UP
hnames[which(hnames == "HALLMARK_UV_RESPONSE_DN")] <- "HALLMARK_UV_RESPONSE"
hprobes[[which(hnames == "HALLMARK_UV_RESPONSE")]] <- c(hprobes[[which(hnames == "HALLMARK_UV_RESPONSE")]],
                                                        hprobes[[which(hnames == "HALLMARK_UV_RESPONSE_UP")]])
hprobes <- hprobes[-which(hnames == "HALLMARK_UV_RESPONSE_UP")]
hnames <- hnames[-which(hnames == "HALLMARK_UV_RESPONSE_UP")]

length(hprobes) # 46
length(hnames) # 46

# Remove "HALLMARK_" from the hnames
hnames_new <- sapply(hnames, function(x){ strsplit(x, "HALLMARK_")[[1]][2] }, USE.NAMES = FALSE)
hnames <- hnames_new

#### Explore which genes belong to which Hallmark ----

## function to ask if a gene symbol is within a particular gene set ##
## returns 1 if in gene set, 0 otherwise ##
pyn <- function(gsym, glist)
{  out <- as.numeric(gsym %in% glist)
return(out) }

## Loop though our diff genes, which Hallmark do they belong to?
hallmark_out <- NULL
for(i in 1:nrow(diff_genes)){
  gname <- diff_genes$GeneName[i]
  ind <- sapply(hprobes, pyn, gsym=gname)
  names(ind) <- hnames
  hallmark_out <- rbind(hallmark_out, ind)	
}
rownames(hallmark_out) <- diff_genes$ProbeID

dim(hallmark_out) # 793 x 46 = genes x broader hallmarks
any_hallmark <- apply(hallmark_out, 1, function(x) { any(x == 1) })
table(any_hallmark)
# FALSE  TRUE 
# 465   328

## Add Hallmarks to table diff_genes --- 
# identical(rownames(hallmark_out), diff_genes$ProbeID) # TRUE
diff_genes$Hallmark <- rep(NA, nrow(diff_genes))
for(i in 1:nrow(hallmark_out)){
  l <- length(which(hallmark_out[i,] == 1))
  if(l > 0){
    hallmarks_temp <- colnames(hallmark_out)[which(hallmark_out[i,] == 1)]
    diff_genes$Hallmark[i] <- paste0(hallmarks_temp[1], ";")
    if(l > 1){
      for(j in 2:length(hallmarks_temp)){
        diff_genes$Hallmark[i] <- c(paste0(diff_genes$Hallmark[i], hallmarks_temp[j], ";"))
      }
    }
  }
}


#### START GO BP ----------
# We also use gene sets from the GO BP dataset, which corresponds to similar biological processes as the Hallmarks
# Find which genes belongs to which GO gene groups

all(c5names_use %in% hnames) # TRUE
all(hnames %in% c5names_use) # TRUE

c5names <- c5names_use
c5probes <- c5probes_use
length(c5names) # 46
length(c5probes) # 46

# Now, do the same thing as with Hallmarks: 

## Loop though our diff genes, in which GO groups are they in?
GO_out <- NULL
for(i in 1:nrow(diff_genes)){
  gname <- diff_genes$GeneName[i]
  ind <- sapply(c5probes, pyn, gsym=gname)
  names(ind) <- c5names
  GO_out <- rbind(GO_out, ind)	
}
rownames(GO_out) <- diff_genes$ProbeID

dim(GO_out) # 793 x 46 = genes x GO
any_GO <- apply(GO_out, 1, function(x) { any(x == 1) })
table(any_GO)
# FALSE  TRUE 
# 505   288

# any GO or HALLMARK?
out <- cbind(any_hallmark, any_GO)
any_hallmark_or_GO <- apply(out, 1, function(x) { any(x == TRUE) })
table(any_hallmark_or_GO)
# 371   422

## Add GO groups to table diff_genes --- 
# identical(rownames(GO_out), diff_genes$ProbeID) # TRUE
diff_genes$GO_BP <- rep(NA, nrow(diff_genes))
for(i in 1:nrow(GO_out)){
  l <- length(which(GO_out[i,] == 1))
  if(l > 0){
    GO_temp <- colnames(GO_out)[which(GO_out[i,] == 1)]
    diff_genes$GO_BP[i] <- paste0(GO_temp[1], ";")
    if(l > 1){
      for(j in 2:length(GO_temp)){
        diff_genes$GO_BP[i] <- c(paste0(diff_genes$GO_BP[i], GO_temp[j], ";"))
      }
    }
  }
}


##### Add broader groups ##### 
# We merged some of the groups together to broader groups
# They are related to similar functions, and will be together in the Heatmap.

CELL_CYCLE <- c("E2F_TARGETS",
                "G2M_CHECKPOINT",
                "MITOTIC_SPINDLE")
EPITHELIAL_STRUCTURE <- c("APICAL_JUNCTION",
                          "APICAL_SURFACE")
IMMUNE <- c("ALLOGRAFT_REJECTION",
            "COMPLEMENT",
            "INTERFERON_GAMMA_RESPONSE",
            "INTERFERON_ALPHA_RESPONSE",
            "IL2_STAT5_SIGNALING",
            "IL6_JAK_STAT3_SIGNALING",
            "TNFA_SIGNALING_VIA_NFKB",
            "TGF_BETA_SIGNALING",
            "INFLAMMATORY_RESPONSE")
METABOLIC <- c("BILE_ACID_METABOLISM",
               "CHOLESTEROL_HOMEOSTASIS",
               "FATTY_ACID_METABOLISM",
               "GLYCOLYSIS",
               "HEME_METABOLISM",
               "OXIDATIVE_PHOSPHORYLATION",
               "PEROXISOME",
               "XENOBIOTIC_METABOLISM")
PI3K_AKT_MTOR <- c("MTORC1_SIGNALING",
                   "PI3K_AKT_MTOR_SIGNALING")

# Add broader Hallmarks to diff_genes
diff_genes$Hallmark_broader <- rep("", nrow(diff_genes))

for(i in 1:nrow(diff_genes)){
  hallmarks_temp <- unlist(strsplit(diff_genes$Hallmark[i], ";"))
  if(any(hallmarks_temp %in% CELL_CYCLE)){ diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], "CELL_CYCLE;")  }
  if(any(hallmarks_temp %in% EPITHELIAL_STRUCTURE)){ diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], "EPITHELIAL_STRUCTURE;")  }
  if(any(hallmarks_temp %in% IMMUNE)){ diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], "IMMUNE;")  }
  if(any(hallmarks_temp %in% METABOLIC)){ diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], "METABOLIC;")  }
  if(any(hallmarks_temp %in% PI3K_AKT_MTOR)){ diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], "PI3K_AKT_MTOR;")  }
  # Remove the ones in the broader groups
  r <- which(hallmarks_temp %in% CELL_CYCLE | hallmarks_temp %in% EPITHELIAL_STRUCTURE
             | hallmarks_temp %in% IMMUNE | hallmarks_temp %in% METABOLIC | hallmarks_temp %in% PI3K_AKT_MTOR | is.na(hallmarks_temp))
  if(length(r) > 0){ hallmarks_temp <- hallmarks_temp[-r] }
  # And add the rest, that are not in the broader groups
  if(length(hallmarks_temp) > 0){
    for(j in 1:length(hallmarks_temp)) {
      diff_genes$Hallmark_broader[i] <- paste0(diff_genes$Hallmark_broader[i], hallmarks_temp[j], ";")
    }
  }
}
diff_genes$Hallmark_broader[which(diff_genes$Hallmark_broader == "")] <- NA

# Add broader GO to diff_genes
diff_genes$GO_broader <- rep("", nrow(diff_genes))

for(i in 1:nrow(diff_genes)){
  GO_temp <- unlist(strsplit(diff_genes$GO_BP[i], ";"))
  if(any(GO_temp %in% CELL_CYCLE)){ diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], "CELL_CYCLE;")  }
  if(any(GO_temp %in% EPITHELIAL_STRUCTURE)){ diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], "EPITHELIAL_STRUCTURE;")  }
  if(any(GO_temp %in% IMMUNE)){ diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], "IMMUNE;")  }
  if(any(GO_temp %in% METABOLIC)){ diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], "METABOLIC;")  }
  if(any(GO_temp %in% PI3K_AKT_MTOR)){ diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], "PI3K_AKT_MTOR;")  }
  # Remove the ones in the broader groups
  r <- which(GO_temp %in% CELL_CYCLE | GO_temp %in% EPITHELIAL_STRUCTURE
             | GO_temp %in% IMMUNE | GO_temp %in% METABOLIC | GO_temp %in% PI3K_AKT_MTOR | is.na(GO_temp))
  if(length(r) > 0){ GO_temp <- GO_temp[-r] }
  # And add the rest, that are not in the broader groups
  if(length(GO_temp) > 0){
    for(j in 1:length(GO_temp)) {
       diff_genes$GO_broader[i] <- paste0(diff_genes$GO_broader[i], GO_temp[j], ";")
    }
  }
}
diff_genes$GO_broader[which(diff_genes$GO_broader == "")] <- NA


##### Now categorize all genes, depending on their Hallmarks broader groups / GO BP broader groups #####

#### We want to add Gene groups to the genes ----
diff_genes$Group <- rep(NA, nrow(diff_genes))

#### THE RULES ####

# 1) Histone or Homeobox ----

## Add HISTONE
# diff_genes$GeneName[grep("HIST", diff_genes$GeneName)]
diff_genes$Group[grep("HIST", diff_genes$GeneName)] <- "HISTONE"

## Add HOMEOBOX
homeobox_genes <- c(grep("HOX", diff_genes$GeneName), 
                    grep("DLX", diff_genes$GeneName), 
                    grep("IRX", diff_genes$GeneName), 
                    grep("NKX", diff_genes$GeneName), 
                    grep("SIX", diff_genes$GeneName))
# diff_genes[homeobox_genes,]
diff_genes$Group[homeobox_genes] <- "HOMEOBOX"

# 2) If one hallmark_broader group, use this one ----

# Find the genes with only one Hallmark_broader group (and not NA)
multiple_hallmarks <- sapply(diff_genes$Hallmark_broader, function(x){ length(strsplit(x,";")[[1]]) > 1 })
genes_to_remove <- c(which(multiple_hallmarks == TRUE),  which(is.na(diff_genes$Hallmark_broader)))
genes_of_interest <- diff_genes$GeneName[-genes_to_remove]
length(genes_of_interest) # 194
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 5
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 189
# Add Hallmark_broader group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- unlist(strsplit(diff_genes$Hallmark_broader[these], ";"))

length(which(!is.na(diff_genes$Group))) # 247
length(which(is.na(diff_genes$Group)))  # 546

# 3) if no hallmark group, but one GO_broader, use this one ----

# Find the genes with only one GO_broader group, and none Hallmark_broader group
multiple_GO <- sapply(diff_genes$GO_broader, function(x){ length(strsplit(x,";")[[1]]) > 1 })
genes_to_remove <- c(which(multiple_GO == TRUE),  which(is.na(diff_genes$GO_broader)), which(!is.na(diff_genes$Hallmark_broader)))
genes_of_interest <- diff_genes$GeneName[-genes_to_remove]
length(genes_of_interest) # 70
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 13
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 57
# Add GO BP group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- unlist(strsplit(diff_genes$GO_broader[these], ";"))

length(which(!is.na(diff_genes$Group))) # 304
length(which(is.na(diff_genes$Group)))  # 489

# 4) If Hallmark_broader and GO_broader group overlaps, select this. ----

# Find the genes with any similar genes in Hallmark_broader and GO_broader groups
genes_of_interest <- c()
group_for_genes_of_interest <- c()
for(i in 1:nrow(diff_genes)){
  if( !is.na(diff_genes$Hallmark_broader[i]) & !is.na(diff_genes$GO_broader[i]) ){
    # select the groups
    if( multiple_hallmarks[i] == FALSE ){ h <- strsplit(diff_genes$Hallmark_broader[i], ";")[[1]]
    }else{ h <- unlist(strsplit(diff_genes$Hallmark_broader[i], ";")) }
    if( multiple_GO[i] == FALSE ){ g <- strsplit(diff_genes$GO_broader[i], ";")[[1]]
    }else{ g <-  unlist(strsplit(diff_genes$GO_broader[i], ";")) }
    # check if any similar, and only 1 similar!
    if(any(g %in% h) & length(which(g %in% h)) == 1 ) {
      genes_of_interest <- c(genes_of_interest, diff_genes$GeneName[i])
      group_for_genes_of_interest <- c(group_for_genes_of_interest, g[which(g %in% h)])
      }
  }
}
length(genes_of_interest) # 92
length(group_for_genes_of_interest) # 92
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 40
group_for_genes_of_interest <- group_for_genes_of_interest[-which(genes_of_interest %in% in_group)]
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 52
length(group_for_genes_of_interest) # 52
# Add group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- group_for_genes_of_interest

length(which(!is.na(diff_genes$Group))) # 356
length(which(is.na(diff_genes$Group)))  # 437

## Look up which groups are the really small (< 3 in each)
hnames_wo_broader <- hnames[-which(hnames %in% CELL_CYCLE | hnames %in% EPITHELIAL_STRUCTURE
                                   | hnames %in% IMMUNE | hnames %in% METABOLIC | hnames %in% PI3K_AKT_MTOR)]
t <- factor(diff_genes$Group, levels=c(hnames_wo_broader, "CELL_CYCLE", "EPITHELIAL_STRUCTURE", "IMMUNE",
                                       "METABOLIC", "PI3K_AKT_MTOR", "HISTONE", "HOMEOBOX"))
t <- table(t)
t <- t[order(t, decreasing=TRUE)]
t <- as.data.frame(t)

really_small_groups <- as.character(t$t[which(t$Freq < 3)])

#### Now repeat, without the small groups ####

# Add Hallmark groups (without small groups) to diff_genes
diff_genes$Hallmark_wo_small <- rep("", nrow(diff_genes))

for(i in 1:nrow(diff_genes)){
  hallmarks_temp <- unlist(strsplit(diff_genes$Hallmark_broader[i], ";"))
  
  if(any(hallmarks_temp %in% really_small_groups)){
    hallmarks_temp <- hallmarks_temp[-which(hallmarks_temp %in% really_small_groups)]
  }
  if(any(is.na(hallmarks_temp))) { hallmarks_temp <- hallmarks_temp[-which(is.na(hallmarks_temp))] }
  if(length(hallmarks_temp) > 0){
    for(j in 1:length(hallmarks_temp)) {
      diff_genes$Hallmark_wo_small[i] <- paste0(diff_genes$Hallmark_wo_small[i], hallmarks_temp[j], ";")
    }
  }
}
diff_genes$Hallmark_wo_small[which(diff_genes$Hallmark_wo_small == "")] <- NA

# Add GO groups (without small groups) to diff_genes
diff_genes$GO_wo_small <- rep("", nrow(diff_genes))

for(i in 1:nrow(diff_genes)){
  GO_temp <- unlist(strsplit(diff_genes$GO_broader[i], ";"))
  
  if(any(GO_temp %in% really_small_groups)){
    GO_temp <- GO_temp[-which(GO_temp %in% really_small_groups)]
  }
  if(any(is.na(GO_temp))) { GO_temp <- GO_temp[-which(is.na(GO_temp))] }
  if(length(GO_temp) > 0){
    for(j in 1:length(GO_temp)) {
      diff_genes$GO_wo_small[i] <- paste0(diff_genes$GO_wo_small[i], GO_temp[j], ";")
    }
  }
}
diff_genes$GO_wo_small[which(diff_genes$GO_wo_small == "")] <- NA

# 5) If one hallmark group - use this one ----

# Find the genes with only one Hallmark_wo_small group
multiple_hallmarks <- sapply(diff_genes$Hallmark_wo_small, function(x){ length(strsplit(x,";")[[1]]) > 1 })
genes_to_remove <- c(which(multiple_hallmarks == TRUE),  which(is.na(diff_genes$Hallmark_wo_small)))
genes_of_interest <- diff_genes$GeneName[-genes_to_remove]
length(genes_of_interest) # 213
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 192
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 21
# Add Hallmark_wo_small group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- unlist(strsplit(diff_genes$Hallmark_wo_small[these], ";"))

length(which(!is.na(diff_genes$Group))) # 377
length(which(is.na(diff_genes$Group)))  # 416

# 6) if no Hallmark_wo_small group, but one GO_wo_small - use this one ----

# Find the genes with only one GO_wo_small group, and none Hallmark_wo_small group
multiple_GO <- sapply(diff_genes$GO_wo_small, function(x){ length(strsplit(x,";")[[1]]) > 1 })
genes_to_remove <- c(which(multiple_GO == TRUE),  which(is.na(diff_genes$GO_wo_small)), which(!is.na(diff_genes$Hallmark_wo_small)))
genes_of_interest <- diff_genes$GeneName[-genes_to_remove]
length(genes_of_interest) # 77
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 72
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 5
# Add GO_wo_small group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- unlist(strsplit(diff_genes$GO_wo_small[these],";"))

length(which(!is.na(diff_genes$Group))) # 382
length(which(is.na(diff_genes$Group)))  # 411

# 7) if ESTROGEN_RESPONSE among hallmark_broader or GO_broader groups - select ESTROGEN_RESPONSE ----
ind <- c(grep("ESTROGEN_RESPONSE", diff_genes$GO_BP),  grep("ESTROGEN_RESPONSE", diff_genes$Hallmark))
ind <- unique(ind[order(ind)])
genes_of_interest <-  diff_genes$GeneName[ind]
length(genes_of_interest) # 50
# Remove the ones that already are in a group:
in_group <- genes_of_interest[which(!is.na(diff_genes$Group[which(diff_genes$GeneName %in% genes_of_interest)]))]
length(in_group) # 39
genes_of_interest <- genes_of_interest[-which(genes_of_interest %in% in_group)]
length(genes_of_interest) # 11
# Add group!
these <- which(diff_genes$GeneName %in% genes_of_interest)
diff_genes$Group[these] <- "ESTROGEN_RESPONSE"

length(which(!is.na(diff_genes$Group))) # 393
length(which(is.na(diff_genes$Group)))  # 400

both_NA <- which(is.na(diff_genes$Hallmark) & is.na(diff_genes$GO_broader) & is.na(diff_genes$Group) )
length(both_NA) # 332

test <- diff_genes[-both_NA,]
test2 <- test[which(is.na(test$Group)), ]
nrow(test2) # 68 left to categorize, this is very much OK!

hnames_wo_broader <- hnames[-which(hnames %in% CELL_CYCLE | hnames %in% EPITHELIAL_STRUCTURE
                                   | hnames %in% IMMUNE | hnames %in% METABOLIC | hnames %in% PI3K_AKT_MTOR)]
t <- factor(diff_genes$Group, levels=c(hnames_wo_broader, "CELL_CYCLE", "EPITHELIAL_STRUCTURE", "IMMUNE",
                                        "METABOLIC", "PI3K_AKT_MTOR", "HISTONE", "HOMEOBOX"))
t <- table(t)
t <- t[order(t, decreasing=TRUE)]
t <- as.data.frame(t)
t
#                                    t Freq
# 1                         CELL_CYCLE   68
# 2                             IMMUNE   60
# 3                          METABOLIC   51
# 4                            HISTONE   45
# 5                  ESTROGEN_RESPONSE   30
# 6                          APOPTOSIS   20
# 7                         DNA_REPAIR   17
# 8                      PI3K_AKT_MTOR   14
# 9               EPITHELIAL_STRUCTURE   13
# 10                          HOMEOBOX   13
# 11                       MYC_TARGETS    9
# 12 EPITHELIAL_MESENCHYMAL_TRANSITION    8
# 13                 PROTEIN_SECRETION    7
# 14   REACTIVE_OXIGEN_SPECIES_PATHWAY    7
# 15                    KRAS_SIGNALING    7
# 16                       P53_PATHWAY    6
# 17        WNT_BETA_CATENIN_SIGNALING    3
# 18                           HYPOXIA    2
# 19                      ADIPOGENESIS    2
# 20         UNFOLDED_PROTEIN_RESPONSE    2
# 21                       UV_RESPONSE    2
# 22                   SPERMATOGENESIS    2
# 23                   NOTCH_SIGNALING    1
# 24                 ANDROGEN_RESPONSE    1
# 25                        MYOGENESIS    1
# 26                      ANGIOGENESIS    1
# 27                       COAGULATION    1
# 28                HEDGEHOG_SIGNALING    0
# 29               PANCREAS_BETA_CELLS    0


### Save

diff_genes$Ultralow_expression <- ifelse(diff_genes$tstat < 0, "Down-regulated", "Up-regulated")

diff_genes_for_paper <- diff_genes[,c("ProbeID", "GeneName", "tstat", "Ultralow_expression", "pvalue", "FDR", "sens", "Group")]
# Fix t-stat
diff_genes_for_paper$tstat <- formatC(diff_genes_for_paper$tstat, 2, format="f")
# Fix p-value
diff_genes_for_paper$pvalue <- formatC(diff_genes_for_paper$pvalue, 2, format="E")
diff_genes_for_paper$pvalue[which(diff_genes_for_paper$pvalue == "0.00E+00")] <- "<1.00E-07"
# Fix FDR
diff_genes_for_paper$FDR <- formatC(diff_genes_for_paper$FDR, 2, format="E")
diff_genes_for_paper$FDR[which(diff_genes_for_paper$FDR == "0.00E+00")] <- "<1.00E-07"
# Fix sensitivity
diff_genes_for_paper$sens <- formatC(diff_genes_for_paper$sens, 3, format="f")
diff_genes_for_paper$sens[which(diff_genes_for_paper$sens == "0.000")] <- "<0.001"
# Fix group
diff_genes_for_paper$Group[which(is.na(diff_genes_for_paper$Group))] <- ""
# New colnames
colnames(diff_genes_for_paper) <- c("ProbeID", "Gene",  "t-statistics", "Ultralow_expression", "p-value", "FDR", "Sensitivity", "Group")

# write.table(diff_genes, "eTable2_full.txt", sep="\t", row.names=F, quote=F)
# write.table(diff_genes_for_paper, "eTable2_for_paper.txt", sep="\t", row.names=F, quote=F)



