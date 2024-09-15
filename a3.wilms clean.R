
################################################################################
Packages 

#if (!require("BiocManager"))
# install.packages("BiocManager")
#BiocManager::install("maEndToEnd")

suppressPackageStartupMessages({library("maEndToEnd")})

#install.packages("Biostrings")  
#install.packages("GenomeInfoDb")

#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(pd.hg.u133a)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)

################################################################################

#Downloaded data via FTP

#import the sdfr file
sdrf_location = file.path(
  "C:/Users/hazel/OneDrive/Documents/UoL/Analytical Skills in Precision Medicine/assessment_2_R/E-GEOD-31403.sdrf.txt")
SDRF = read.delim(sdrf_location)

#sample names (Array.Data.File) used as rownames 
rownames(SDRF) = SDRF$Array.Data.File

#turn the SDRF table into an AnnotatedDataFrame, will need later to create 
#an ExpressionSet for the data
SDRF = AnnotatedDataFrame(SDRF)

#import CEL files: checks if CEL files in the order that 
#corresponds to the SDRF table 
raw_data = oligo::read.celfiles(filenames = file.path("C:/Users/hazel/OneDrive/Documents/UoL/Analytical Skills in Precision Medicine/assessment_2_R",
                                                      SDRF$Array.Data.File),
                                verbose = FALSE, phenoData = SDRF)

# to check the objects are valid
stopifnot(validObject(raw_data))

#pData used to directly access the phenoData in the expressionSet raw_data
head(Biobase::pData(raw_data))

# subselect the wanted columns
Biobase::pData(raw_data) = Biobase::pData(raw_data)[,c("Hybridization.Name",
                                                       "Characteristics..Subset.")]

################################################################################

#Quality control/assessment of the raw data

#The expression intensity values are in the assayData sub-object “exprs”
#can be accessed by the exprs(raw_data) function. 
#The rows represent the microarray probes, 
#i.e.  the single DNA locations on the chip, while the 
#columns represent one microarray

hhh = Biobase::exprs(raw_data)[1:5, 1:5]
nrows(hhh)

length(raw_data@assayData$exprs)
length(raw_data@featureData@data)
raw_data@featureData@data
raw_data@assayData$exprs

# take the log2 of Biobase::exprs(raw_data), as expression data is commonly 
#analysed on a logarithmic scale

raw_exp = log2(Biobase::exprs(raw_data))

#perform a PCA
raw_PCA = prcomp(t(raw_exp), scale. = FALSE)

percenVar = round(100*raw_PCA$sdev^2/sum(raw_PCA$sdev^2),1)

sd_ratio = sqrt(percenVar[2]/percenVar[1])

#PCA to show if stage of CA or if tumour is from relapse is causing variance
#neither are
dataGG = data.frame(PC1 = raw_PCA$x[,1], PC2 = raw_PCA$x[,2],
                    Subset = pData(raw_data)$Characteristics..Subset.)

dataGG$Subset = as.factor(dataGG$Subset)

ggplot(dataGG,aes(PC1, PC2)) +
  geom_point(aes(colour = Subset))+
  ggtitle("PCA plot of the log-transformed pre-processed array data") +
  xlab(paste0("PC1, VarExp: ", percenVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percenVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("dodgerblue4","darkgreen","red3","orange","purple3"))


oligo::boxplot(raw_data, target = "core", las = 2, xaxt = "n",
               main = "Boxplot of log2-intensitites of the pre-processed array data",
               xlab = "Arrays",
               ylab = "Intensity")
#intensity distributions of the individual arrays are quite different, 
#indicating the need for an appropriate normalization


#The package produces an html report, containing the quality control plots 
#together with a description of their aims and an identification of possible 
#outliers. 

#arrayQualityMetrics(expressionset = raw_data,
#                    outdir = tempdir(),
#                    force = TRUE, do.logtransform = TRUE,
#                   intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))


#Background adjustment
#After the initial import and quality assessment, the next step in processing of
#microarray data is background adjustment. This is essential because a proportion 
#of the measured probe intensities are due to non-specific hybridization and the
#noise in the optical detection system. Therefore, observed intensities need to
#be adjusted to give accurate measurements of specific hybridization.


#Across-array normalization (calibration)
#Normalization across arrays is needed in order to be able to compare 
#measurements from different array hybridizations due to many obscuring sources 
#of variation. These include different efficiencies of reverse transcription, 
#labeling or hybridization reactions, physical problems with the arrays, 
#reagent batch effects, and laboratory conditions.

#Summaraization 
#After normalization, summarization is necessary to be done because on the 
#Affymetrix platform, transcripts are represented by multiple probes, that is 
#multiple locations on the array. For each gene, the background-adjusted and 
#normalized intensities of all probes need to be summarized into one quantity 
#that estimates an amount proportional to the amount of RNA transcript.

#One-step preprocessing 
#The package oligo allows us to perform background correction, normalization and
#summarization in one single step using a deconvolution method for 
#background correction, quantile normalization and the 
#RMA (robust multichip average) algorithm for summarization.

#This series of steps as a whole is commonly referred to as RMA algorithm, 
#although strictly speaking RMA is merely a summarization method 

################################################################################

#RMA calibration of the data 

#full RMA algorithm to our data in order to background-correct, 
#normalize and summarize

#library(affy)

data_rma_norm = oligo::rma(raw_data)

#Quality assessment of the calibrated data

data_rma_norm_exp = log2(Biobase::exprs(data_rma_norm))

#perform a PCA
data_rma_norm_PCA = prcomp(t(data_rma_norm_exp), scale. = FALSE)

percentVar = round(100*data_rma_norm_PCA$sdev^2/sum(data_rma_norm_PCA$sdev^2),1)
sd_ratio = sqrt(percentVar[2] / percentVar[1])

PCA_data_rma = data.frame(PC1 = data_rma_norm_PCA$x[,1], PC2 = data_rma_norm_PCA$x[,2],
                          Subset =data_rma_norm$Characteristics..Subset.)

# Replace values in 'column_name'
#Need to change individually
PCA_data_rma$Subset = ifelse(PCA_data_rma$Subset 
                             == "1", "1", PCA_data_rma$Subset)

ggplot(PCA_data_rma, aes(PC1, PC2)) +
  geom_point(aes(colour = Subset)) +
  ggtitle("PCA plot of the normalised array data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("dodgerblue4","darkgreen","red3","orange",
                                "purple3"))

#should show more variance in PC1 than PC2, compared to raw PCA


oligo::boxplot(data_rma_norm, target = "core", las = 2, xaxt = "n",
               main = "Boxplot of log2-intensitites for the normalised array data",
               xlab = "Arrays",
               ylab = "Intensity")

################################################################################

#Filtering based on intensity 

#We now filter out lowly expressed genes.
#For intensity-based filtering, calculate the row-wise medians from the 
#expression data, as they represent the transcript medians

data_medians = rowMedians(Biobase::exprs(data_rma_norm))

hist_res = hist(data_medians, 100, col = "palevioletred3", freq = FALSE, 
                main = "Histogram of the median intensities", 
                border = "antiquewhite4",
                xlab = "Median intensities")


#In the histogram of the gene-wise medians, can see an enrichment of 
#low medians on the left hand side. 
#These represent the genes to filter. 
#In order to infer a cutoff from the data, inspect the histogram: 
#visually set a cutoff line man_threshold to the left of the histogram peak in 
#order not to exclude too many genes.
man_threshold = 6

hist_res = hist(data_medians, 100, col = "lightgoldenrod2", freq = FALSE, 
                main = "Histogram of the median intensities", 
                border = "antiquewhite4",
                xlab = "Median intensities")

abline(v = man_threshold, col = "royalblue2", lwd = 3)


#Transcripts that do not have intensities larger than the threshold in at least 
#as many arrays as the smallest experimental group are excluded.

#first have to get a list with the number of samples
no_of_samples = 
  table(paste0(pData(data_rma_norm)$Characteristics..Subset.
                              ,"_",pData(data_rma_norm)$FactorValue..RELAPSE.
  ))
no_of_samples 

#filter out all transcripts that do not have intensities greater than the 
#threshold in at least as many arrays as the smallest experimental group (2) 
#(samples_cutoff).

#A function idx_man_threshold is then applied to each row, 
#i.e.  to each transcript across all arrays. 
#It evaluates whether the number of arrays where the median intensity passes 
#the threshold (sum(x > man_threshold)) is greater than the samples_cutoff and 
#returns TRUE or FALSE for each row, i.e. each transcript.

samples_cutoff = min(no_of_samples)
samples_cutoff


idx_man_threshold = apply(Biobase::exprs(data_rma_norm), 1,
                          function(x){
                            sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

data_manfiltered = subset(data_rma_norm, idx_man_threshold)

################################################################################

#Annotation of the transcript clusters

#library(AnnotationDbi)

#BiocManager::install("hgu133a.db")

library(hgu133a.db)

#annotation information to the transcript cluster identifiers stored in the 
#featureData of our ExpressionSet
probe_ids = featureNames(data_manfiltered)

data_anno = AnnotationDbi::select(hgu133a.db,
                                  keys = probe_ids,
                                  columns = c("SYMBOL", "GENENAME"),
                                  keytype = "PROBEID")
data_anno

#filtered out the probes that do not map to a gene
data_anno = subset(data_anno, !is.na(SYMBOL))

# compute a summary table to see how many transcript-cluster identifiers 
#had mapped to multiple gene symbols

dat_anno_grouped = group_by(data_anno, PROBEID)
data_anno_summarised = 
  dplyr::summarize(dat_anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(data_anno_summarised)

data_anno_filtered = filter(data_anno_summarised, no_of_matches > 1)

head(data_anno_filtered)

probe_stats = data_anno_filtered 

nrow(probe_stats)
#575 probes mapped to more than 1 transcript-cluster identifier

dim(probe_stats)


#exclude probe ids that match to more than one gene
ids_to_exlude = (featureNames(data_manfiltered) %in% probe_stats$PROBEID)

table(ids_to_exlude)

data_final = subset(data_manfiltered, !ids_to_exlude)

validObject(data_final)

#exclude probe ids from the feature data
#fData enables the to access the feature data of an expression set. 
#Until now, no feature data whatsoever is stored in the fData(data_final). 
#Only the row names are the row names of the assay data by default, 
#which are the PROBEIDs of the transcripts.
#generate a column PROBEID in fData(data_final) and assign the row names of 
#fData(data_final)

#Therefore, we generate a column PROBEID in fData(data_final) and 
#assign the row names of fData(data_final) to it
fData(data_final)$PROBEID = rownames(fData(data_final))

#left-join fData(data_final)with data_anno, which already contains the 
#columns “SYMBOL” and “GENENAME”. 
#A left-join keeps the rows and columns of the first argument and 
#adds the corresponding column entries of the second argument
fData(data_final) = left_join(fData(data_final), data_anno)

#Restoring rnames after left-join
rownames(fData(data_final)) = fData(data_final)$PROBEID 

validObject(data_final)

################################################################################

#Graphs of the filtered data before any stat testing

#Nice PCA
percentVar.final = round(100*PCA_final$sdev^2/sum(PCA_final$sdev^2),1)
sd_ratio.final = sqrt(percentVar.final[2] / percentVar.final[1])

PCA_final.df = data.frame(PC1 = PCA_final$x[,1], PC2 = PCA_final$x[,2],
                          Subset =data_final$Characteristics..Subset.)

# Replace values in 'column_name'
#Need to change individually
PCA_final.df$Subset = ifelse(PCA_final.df$Subset 
                             == "5", "5", PCA_final.df$Subset)

final_PCA_plot = ggplot2::ggplot(PCA_final.df, aes(PC1, PC2)) +
  geom_point(aes(colour = Subset)) +
  ggtitle("PCA plot of the filtered data") +
  xlab(paste0("PC1, VarExp: ", percentVar.final[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar.final[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("dodgerblue4","darkgreen","red3","orange",
                                "purple3"))

final_PCA_plot


oligo::boxplot(data_final, target = "core", las = 2, xaxt = "n",
               main = "Boxplot of log2-intensitites for the filtered, normalised data",
               xlab = "Arrays",
               ylab = "Intensity")

################################################################################

#Linear models S5 vs S1-4 and ebayes analysis

head(pData(data_final))

#To find differential expression between subset 5 and 1-4
#Experimental variables: subset 5 and subset 1-4
individual = 
  as.character(Biobase::pData(data_final)$Hybridization.Name)

Subset = ifelse(str_detect(Biobase::pData(data_final)$Characteristics..Subset., 
                           
                           "5"), "S5", "S1_and_S4")

#Creating design matrix 
#~0 means no intercept term in the model 
data_design_S5 = model.matrix(~ 0 + Subset)

rownames(data_design_S5) = individual  

#Can do single gene analysis

#contrast and hypothesis tests

#fit a linear model to the data and apply the contrasts.fit() function to it in 
#order to find genes with significant differential expression 
#between S5 and S1-4

data_fit_S5 = lmFit(data_final,design = data_design_S5)

??makeContrasts

contrast_matrix_S5 = makeContrasts(SubsetS1_and_S4-SubsetS5,
                                   levels = data_design_S5)

contrast_matrix_S5


contrasts_result = contrasts.fit(data_fit_S5, contrasts = contrast_matrix_S5)
contrasts_result

#applied the empirical Bayes variance moderation method to the model via 
#the eBayes() function, which computes moderated t-statistics. 
#Using a combination of the per-gene-variance and a prior variance we can improve
#the variance estimate, hence the term “moderation”. “Empirical Bayes” means 
#that the prior is estimated from the data.
#The result of the eBayes() step is that the individual variances are shrunken 
#towards the prior value.

data_model_S5 = eBayes(contrasts_result)

#Extracting results

#to extract the results from the t test 
table_S5 = topTable(data_model_S5, number = Inf)
head(table_S5)
table_S5[1:10,]

dim(table_S5)

summary(isUnique(table_S5$SYMBOL))

#ordering by FDR
table_S5 = table_S5[order(table_S5$adj.P.Val), ]
hh = table_S5[order(table_S5$logFC, decreasing = TRUE), ]
head(hh)

# need to remove repeating genes as volcano shows them and shows in top 10 lists
table_S5_unique = subset(table_S5, !duplicated(SYMBOL))

table_S5_unique = subset(table_S5_unique, !is.na(SYMBOL))

validObject(table_S5_unique)

dim(table_S5_unique)

summary(isUnique(table_S5_unique$SYMBOL))


log.5 = subset(table_S5_unique, table_S5_unique$logFC >= 1  | table_S5_unique$logFC <= -1 )

dim(log.5)


#As a diagnostic check, also plot the p-value histogram: 
#expect a uniform distribution for the p-values that correspond to true null
#hypotheses, while a peak near zero shows an enrichment for low p-values 
#corresponding to differentially expressed (DE) genes.

hist(table_S5_unique$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "subset 5 vs Subset 1 to 4 - Wilm's Tumour", xlab = "p-values")

#total number of differentially expressed genes
total_genenumber_S5 = length(subset(table_S5_unique, adj.P.Val <= 0.05)$SYMBOL)

total_genenumber_S5

#3864
#What was the total number before?
#9944 dim(table_S5_unique)

#Visualising DE analysis results 

#BiocManager::install('EnhancedVolcano')
#library(ggrepel)
#library(EnhancedVolcano)

volcano_plot_result = EnhancedVolcano::EnhancedVolcano(table_S5_unique,
                lab = table_S5_unique$SYMBOL,
                title = "Volcano plot of the differentially expressed genes between 
                Subset 5 and Subsets 1-4.",
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 13e-10,
                FCcutoff = 1,
                pointSize = 2
)

volcano_plot_result

#1.30103 = -log10 of 0.05
-log10(0.05)

top.S5 = table_S5_unique[order(table_S5_unique$adj.P.Val),]
#writeLines(top.S5[1:100,]$SYMBOL, "S5_top_100")

top.S5[1:10,]

top.S5[top.S5$SYMBOL =="STAT3",]



log.5

log_top.S5 = log.5[order(log.5$adj.P.Val),]
log_top.S5[1:10,]

################################################################################

#T statistics 

#T test S5
dat5 = Biobase::exprs(data_final)

t2.5 = vector()

dat.prior.5 = log2(Biobase::exprs(raw_data))

pval.t2.5 = vector()

group = SDRF@data$Characteristics..Subset.


for(j in 1:nrow(dat5)){
  temp = dat5[j,]
  res=t.test(temp[group!=5], temp[group==5], var.equal=T)
  t2.5[j] = res$stat
  pval.t2.5[j]=res$p.val
}

#str(res)
#res$stat
#res$p.val

adj.pval.t2.5 = p.adjust(pval.t2.5, "BH")

result.table.5 = data.frame(ID=rownames(dat5), t.stat=t2.5,
                            pvalue=pval.t2.5, fdr.pvalue=adj.pval.t2.5)

result.table.5 = subset(result.table2.5, 
                        !duplicated(result.table.5))
result.table.5$fdr.pvalue

result.table.sorted.5 = result.table.5[order(result.table.5$fdr.pvalue),]

result.table.sorted.5[1:10,]

#BiocManager::install("GEOquery")
library(GEOquery)

gse = GEOquery::getGEO("GSE31403", GSEMatrix = TRUE)

feature.data = gse$GSE31403_series_matrix.txt.gz@featureData@data
feature.data = feature.data[,c(1,11)]
feature.data = as.data.frame(feature.data)

top.t.5  = result.table.sorted.5 %>% inner_join(
  ., feature.data, by="ID")

top.t.5 = subset(top.t.5, !duplicated(`Gene Symbol`))

#writeLines(top.t.5$`Gene Symbol`[1:100], "S5_top_100_T")

top.t.5[1:10,]

################################################################################

#GSEA S5 vs S1-4

# Extract all significant results for ORA/GSEA

# Rank results for GSEA; put in decreasing order.
GSEA_data.5 = table_S5_unique %>% filter(!is.na(adj.P.Val)) %>%
  arrange(desc(logFC)) %>% 
  # create a named vector with log fold changes
  dplyr::pull(logFC,name=SYMBOL)

GSEA_data.5
table_S5_unique[1,]

#The gene label permutation step uses random reordering. Results may differ 
#between runs if you do not use a random seed generator (e.g., set.seed()).
set.seed(235)

gsea_go.5 = gseGO(geneList = GSEA_data.5,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  seed=TRUE)

class(gsea_go.5)

head(gsea_go.5.ordered$leading_edge)

head(gsea_go.5)

#graphs 

#enrichplot offers several functions specific to GSEA output, especially in 
#regard to visualizing the running score and preranked list (gseaplot() and 
#gseaplot2()). gseaplot2() allows us to visualize several enriched gene sets 
#at once.
enrichplot::gseaplot2(gsea_go.5, geneSetID = c(1:5))


#We can use an upset plot to look at fold distributions over gene sets. 
#This requries the installation of the ggupset package.
library(ggupset)

enrichplot::upsetplot(gsea_go.5) 


#barplot of NES 
#install.packages("forcats")
library(ggplot2)
library(forcats)

ggplot2::ggplot(gsea_go.5, 
                showCategory=10, 
                ggplot2::aes(NES, fct_reorder(Description, NES),
                             fill=qvalue)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="blue",size=1)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw() + 
  theme(text=element_text(size=8))+
  xlab("Normalised Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA with GO")


#graphs 
gsea_go.5 = pairwise_termsim(gsea_go.5)

emapplot(gsea_go.5, showCategory = 15)
#edge-thickness (that is the line connecting two pathways) is proportional to 
#the number of overlapping genes between two pathways.

dotplot(gsea_go.5, showCategory=12) +
  labs(
    y = "Top 12 most enriched Biological Process GO terms", 
    title = "The Gene Ratio and Gene Counts of the enriched 
       Biological Process GO terms.")

treeplot(gsea_go.5) +
  labs(title = "Functional clustering of the enriched Biological Process GO 
       terms.")


#ordering 
gsea_go.5[1:5,]

gsea_go.5.ordered = gsea_go.5[order(gsea_go.5@result$p.adjust),] 
gsea_go.5.ordered[1:5,]

gene_list = gsea_go.5.ordered$core_enrichment[1:5]

gene_list

# Split the gene symbols by "/"
genes_separated = unlist(strsplit(gene_list, "/"))

s5.gsea.unique = unique(genes_separated)

s5.gsea.unique

writeLines(s5.gsea.unique, "s5_gsea")

gsea_go.5.ordered$Description[1:5]

gsea_go.5.ordered[1:5,]

################################################################################
#enrichR

if(!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}

library(enrichR)

gene_list_enrich = table_S5_unique$SYMBOL

websitelive = getOption("enrichR.live")
if(websitelive) {
  listEnrichrSites()
  setEnrichrSite("EnrichR")
}

if(websitelive)dbs = listEnrichrDbs()
if(websitelive) head(dbs)

dbs$libraryName

enrich_5 = enrichr(gene_list_enrich, databases = "GO_Biological_Process_2023")

class(enrich_5)

enrich_5_df = data.frame(enrich_5$GO_Biological_Process_2023)

#graphs
enrichplot::gseaplot2(enrich_5_df, geneSetID = c(1:2))


#We can use an upset plot to look at fold distributions over gene sets. 
#This requries the installation of the ggupset package.
library(ggupset)

enrichplot::upsetplot(enrich_5_df) 


#barplot of NES 
#install.packages("forcats")
library(ggplot2)
library(forcats)

ggplot2::ggplot(enrich_5_df, 
                showCategory=10, 
                ggplot2::aes(NES, fct_reorder(Description, NES),
                             fill=qvalue)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="blue",size=1)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw() + 
  theme(text=element_text(size=8))+
  xlab("Normalised Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA with GO")


#graphs 
enrich_5_df = pairwise_termsim(enrich_5_df)

emapplot(enrich_5_df, showCategory = 15)
#edge-thickness (that is the line connecting two pathways) is proportional to 
#the number of overlapping genes between two pathways.

dotplot(enrich_5_df, showCategory=12) +
  labs(
    y = "Top 12 most enriched Biological Process GO terms", 
    title = "The Gene Ratio and Gene Counts of the enriched 
       Biological Process GO terms.")

treeplot(enrich_5_df) +
  labs(title = "Functional clustering of the enriched Biological Process GO 
       terms.")

plotEnrich(enrich_5_df, showTerms = 12, orderBy = "Adjusted.P.value")

################################################################################

GSEA_data.log5 = log_top.S5 %>% filter(!is.na(adj.P.Val)) %>%
  arrange(desc(logFC)) %>% 
  # create a named vector with log fold changes
  dplyr::pull(logFC,name=SYMBOL)

GSEA_data.5
table_S5_unique[1,]

#The gene label permutation step uses random reordering. Results may differ 
#between runs if you do not use a random seed generator (e.g., set.seed()).
set.seed(235)

gsea_go.log5 = gseGO(geneList = GSEA_data.log5,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  seed=TRUE)


ggplot2::ggplot(gsea_go.log5, 
                showCategory=10, 
                ggplot2::aes(NES, fct_reorder(Description, NES),
                             fill=qvalue)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="blue",size=1)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw() + 
  theme(text=element_text(size=8))+
  xlab("Normalised Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA with GO")
################################################################################

#Linear models S5 vs S1-4 and ebayes analysis, trying a different direction

head(pData(data_final))


#To find differential expression between subset 5 and 1-4
#Experimental variables: subset 5 and subset 1-4

individual = 
  as.character(Biobase::pData(data_final)$Hybridization.Name)
###Josh note: data_final is my normalised data 

subsets_try = Biobase::pData(data_final)$Characteristics..Subset.
### Josh note: this would be your CAD groups 
subsets_try = factor(subsets_try)
class(subsets_try)
### Josh note: this should say factor 

#Creating design matrix 
#~0 means no intercept term in the model 
data_design_try = model.matrix(~0 + factor(subsets_try))
data_design_try
rownames(data_design_try) = individual  

#contrast and hypothesis tests

#fit a linear model to the data and apply the contrasts.fit() function to it in 
#order to find genes with significant differential expression 
#between S5 and S1-4

colnames(data_design_try) = c("Subset1", "Subset2", "Subset3","Subset4","Subset5")
### Josh note: need to change the column names to match the contrast matrix 

###Josh note: the order should be CAD 1-6 vs CAD 0 (CONTROL), just add the 
###combos you want testing in the contrast matrix code below 
contrast_matrix_try = makeContrasts(Subset1 - Subset5, 
                                    Subset2 - Subset5, 
                                    Subset3 - Subset5,
                                    Subset4 - Subset5,
                                  levels = data_design_try)


data_fit_try = lmFit(data_final,design = data_design_try)

data_fit_try

contrasts_result_try = contrasts.fit(data_fit_try, contrasts = contrast_matrix_try)
contrasts_result_try$cov.coefficients

data_model_try = eBayes(contrasts_result_try)
data_model_try[1:5,]

#Extracting results

#to extract the results from the t test 
table_try = topTable(data_model_try,  number = Inf)
table_10 = topTable(data_model_try,  number = 10) #to get the top 10 sig genes
table_10
head(table_try)
table_try[1:10,]

table_try = subset(table_try,table_try$adj.P.Val <=0.05)

head(table_try)

results = decideTests(data_model_try) # need this for venn diagram to work
vennDiagram(results)

table_try_unique = subset(table_try, !duplicated(SYMBOL))

table_try_unique = subset(table_try_unique, !is.na(SYMBOL))

validObject(table_try_unique)

dim(table_try_unique)

summary(isUnique(table_try_unique$SYMBOL))
table_try_unique[1:50,]

###### not needed
# Extract F-statistics and associated p-values
f_stats = table_try_unique$F
p_values = table_try_unique$P.Value

# Generate histogram for p-values or F-statistics
hist(f_stats, col = "lightblue", main = "Histogram of P-values")
###### not needed 


EnhancedVolcano::EnhancedVolcano(table_try_unique,
                lab = table_try_unique$SYMBOL,
                title = "Volcano plot of the differentially expressed genes between 
                Subset 5 and Subsets 1-4.",
                xlab = "log2 Average Expression",
                x = 'AveExpr',
                y = 'adj.P.Val',
                pCutoff = 10e-32,
                FCcutoff = 2,
                pointSize = 2
)


top.try = table_try[order(table_try_unique$adj.P.Val),]
top.10 = top.try[1:10,]
top.10$P.Value
writeLines(as.character(top.10$AveExpr), "Top_10_try")
top.10
top.try[1:100,]

### GSEA

GSEA_data_try = top.try %>% filter(!is.na(adj.P.Val)) %>%
  arrange(desc(F)) %>% 
  # create a named vector with log fold changes
  dplyr::pull(F,name=SYMBOL)


GSEA_data_try
top.try[1,]

#The gene label permutation step uses random reordering. Results may differ 
#between runs if you do not use a random seed generator (e.g., set.seed()).
set.seed(235)

gsea_go_try = gseGO(geneList = GSEA_data_try,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  seed=TRUE)

head(gsea_go_try)

library(ggplot2)
library(forcats)

ggplot2::ggplot(gsea_go_try, 
                showCategory=10, 
                ggplot2::aes(NES, fct_reorder(Description, NES),
                             fill=qvalue)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="red",size=2)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw() + 
  theme(text=element_text(size=8))+
  xlab("Normalised Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA analysis of the significant DEGs")

gsea_go_try = pairwise_termsim(gsea_go_try)

emapplot(gsea_go_try, showCategory = 15)
#edge-thickness (that is the line connecting two pathways) is proportional to 
#the number of overlapping genes between two pathways.

dotplot(gsea_go_try, showCategory=12) +
  labs(
    y = "Top 12 most enriched Biological Process GO terms", 
    title = "The Gene Ratio and Gene Counts of the enriched 
       Biological Process GO terms.")

treeplot(gsea_go_try) +
  labs(title = "Functional clustering of the enriched Biological Process GO 
       terms.")

cnetplot(gsea_go_try,
         ShowCatergory = 5,
         layout = "kk",
         colourEdge = TRUE,
         cex_label_category = 0.5,
         vortex.label.font = 9,
         color_gene = "darkblue",
         color_category = "goldenrod2"
) +
  labs(
    title = "Relationship between the enriched Cellular Component Go terms and 
    each term's associated differentially expressed genes.")


clusterProfiler::heatplot(gsea_go_try, showCategory = 8)


################################################################################

#T statistic with log fold change that doesn't work
dat5_try = Biobase::exprs(data_final)
group = SDRF@data$Characteristics..Subset.

###start here
# Initialize vectors to store results
t_stats_try = vector()
p_values_try = vector()

##
log_fold_changes_try = matrix(NA, nrow = nrow(dat5_try), ncol = length(unique(group)) - 1)


##miss this 
# Iterate through each subset (except Subset 5)
for (subset_i in setdiff(unique(group), 5)) {
  # Subset-specific data
  subset_data = dat5_try[, group == subset_i]
  
  # Subset 5 data
  subset_5_data = dat5_try[, group == 5]
  
  # Perform t-test
  res = t.test(subset_data, subset_5_data, var.equal = TRUE)
  
  # Store results
  t_stats_try = c(t_stats_try, res$stat)
  p_values_try = c(p_values_try, res$p.val)
  
  # Calculate log-fold change
  log_fold_changes_try[, subset_i] = rowMeans(subset_data) - rowMeans(subset_5_data)
}

log_fold_changes_try

#### use this one instead
for(j in 1:nrow(dat5_try)){
  temp = dat5_try[j,]
  res=t.test(temp[group==c(1,2,3,4)], temp[group==5], var.equal=T)
  t_stats_try[j] = res$stat
  p_values_try[j]=res$p.val
}

### also ignore this 
log_fold_changes_matrix = matrix(NA, nrow = nrow(dat5_try), ncol = 4)

for (j in 1:nrow(dat5_try)) {
  # Subset the data for the current group
  subset_data = dat5_try[, group != 5]
  subset_5_data = dat5_try[, group == 5]
  log_fold_changes_matrix[, j] = log2(subset_data / subset_5_data)
}

length(subset_data)
length(subset_5_data)

# Print the result
print(log_fold_changes_matrix)

### then to here
# Adjust p-values
adj_p_values_try = p.adjust(p_values_try, "BH")

##go pass this 
# Calculate averages of log-fold changes across all subsets
avg_log_fold_changes_try = rowMeans(log_fold_changes_try, na.rm = TRUE)

avg_log_fold_changes_try

###and continue here
# Assuming 'gene_names' is a vector containing gene names
gene_names_try = rownames(dat5_try)

# Create a results data frame
results_table_t_try = data.frame(
  Gene = gene_names_try,
  T_Statistic = t_stats_try,
  P_Value = p_values_try,
  Adjusted_P_Value = adj_p_values_try
)

head(results_table_try)


#more rubbish 
Avg_Log_Fold_Change = avg_log_fold_changes_try


length(gene_names_try)
length(t_stats_try)
length(p_values_try)
length(adj_p_values_try)
length(avg_log_fold_changes_try)


###finially here 
results_table_t_try = subset(results_table_t_try, 
                        !duplicated(results_table_t_try))

results_table_t_try = subset(results_table_t_try,
                             results_table_t_try$Adjusted_P_Value < 0.05)

results_table_t_try_sorted = results_table_t_try[order(
  results_table_t_try$Adjusted_P_Value),]

results_table_t_try_sorted[1:10,]

#BiocManager::install("GEOquery")
library(GEOquery)

gse = GEOquery::getGEO("GSE31403", GSEMatrix = TRUE)

feature.data = gse$GSE31403_series_matrix.txt.gz@featureData@data
feature.data = feature.data[,c(1,11)]
feature.data = as.data.frame(feature.data)

results_table_t_try_sorted$ID = results_table_t_try_sorted$Gene

top.t.try  = results_table_t_try_sorted %>% inner_join(
  ., feature.data, by="ID")

top.t.try = subset(top.t.try, !duplicated(`Gene Symbol`))

top.t.10 = top.t.try[1:10,]
top.t.10

dim(results_table_t_try)

writeLines(as.character(top.t.10$Adjusted_P_Value), "Top_10_T_try")


top.t.try[1:10,]

colnames(results_table_try_sorted)
colnames(feature.data)
results_table_t_try$Gene

EnhancedVolcano::EnhancedVolcano(results_table_t_try,
                                 lab = results_table_t_try$Gene,
                                 title = "Volcano plot of the differentially expressed genes between 
                Subset 5 and Subsets 1-4.",
                xlab = "log2 Average Expression",
                x = 'T_Statistic',
                y = 'Adjusted_P_Value',
                pCutoff = 10e-32,
                FCcutoff = 2,
                pointSize = 2
)

################################################################################

#Session information
gc()

length(getLoadedDLLs())
#125

sessionInfo()

citation("oligo")

#Microarray analysis workflow primarily adapted from Klaus, B. and Reisenauer,S. 2018. PMID: 30101006
