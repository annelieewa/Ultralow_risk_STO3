# Ultralow_risk_STO3
Characteristics of ER-positive ultralow risk tumors in STO-3

<b>Read me file for:
“Clinical and molecular characteristics of ER-positive breast cancer tumors identified as ultralow risk by the 70-gene signature”</b>

Written by Annelie Johansson

<b>DATA FILES</b>

20170404patient_array.RData</br>
Patient clinical data with matched Agilent microarray data for the STO-3 patients

all_modules.RData</br>
Coefficients for the multi-gene modules

h.all.v6.2.symbols.gmt.txt</br>
Hallmarks 

c5.bp.v6.2.symbols.gmt.txt</br>
GO groups


<b>ANALYSIS CODE</b>

<b>Preparation</b>

create_gene_modules_scores.R</br>
Calculates gene module scores, divides the score into tertiles, and categorizes the scores into two groups the worst vs the other two (which depends upon if high or low score is associated with bad or poor prognosis)</br>
Input: 20170404patient_array.RData, all_modules.RData</br>
Output: gene_modules.RData</br>

<b>Comparison of patient/tumor characteristics</b>

fishers_test.R</br>
Performs Fisher’s test for patient and tumor characteristics, including the gene module scores, comparing:</br>
•	Ultralow risk tumors (UL) vs ER-positive tumors (not UL)</br>
•	Ultralow risk tumors vs Luminal A tumors (ER-positive, not UL)</br>
•	Ultralow risk tumors vs Luminal B tumors (ER-positive, not UL)</br>
Input: 20170404patient_array.RData, gene_modules.RData</br>
Output: table1.txt, stable3.txt</br>

<b>Differential gene expression analysis</b>

OCplusDiffAnalysis_probes.R</br>
Differential gene expression analysis (for probes) using function EOC in R package OCplus, a package that uses t-statistics and False Discovery Rates (FDR) as cut offs for determining genes that are expressed at different levels between two groups. We compare ER-positive UL patients vs not UL (total n=538 patients). Using threshold of FDR < 0.001 results in n=793 significant probes in n=706 genes.</br>
Input: 20170404patient_array.RData</br>
Output: stable4.txt</br>


<b>Define and add gene groups (for the heatmap)</b>

explore_GO_BP_classes.R<.br>
In this script we identify GO BP classes with similar functions as the Hallmarks, and save them in groups with same name as the Hallmarks</br>
Input: h.all.v6.2.symbols.gmt.txt, c5.bp.v6.2.symbols.gmt.txt</br>
Output: GO_use.RData (with c5names_use, c5probes_use)</br>

stable4_add_hallmark_GO_groups.R </br>
In this script we identify which genes belongs to which Hallmarks or GO groups, this information is used for the heatmap.</br>
Input: stable4.txt, h.all.v6.2.symbols.gmt.txt, GO_use.RData</br>
Output: stable4_hallmarks_GO.txt</br>


<b>FIGURES CODE</b>

barplot_ultralow.R</br>
Bar plot of gene-modules, Figure 2.</br>
Input: stable3.txt</br>
Output: Figure2.pdf </br>

heatmap_ultralow.R</br>
Code for the heatmap, Figure 3.</br>
Input: 20170404patient_array.RData, stable4_hallmarks_GO.txt</br>
Output: Figure3.pdf</br>


