
BiocManager::install("ggplot2",destdir="F:/software_run/RWorkspace/downloaded_packages")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("amap")
library("ggplot2")
library("BiocParallel")
register(SnowParam(126L,type="MPI"))
setwd("F:/new_disease/3_label_properties/ovary_hyperexpression")
files <- c("ovary.csv","testis.csv","adipose_subcutaneous.csv","adipose_visceral_omentum.csv","adrenal_gland.csv","artery_aorta.csv","artery_coronary.csv","artery_tibial.csv","bladder.csv","brain_amygdala.csv","brain_anterior_cingulate_cortex_ba24.csv","brain_caudate_basal_ganglia.csv","brain_cerebellar_hemisphere.csv","brain_cerebellum.csv","brain_cortex.csv","brain_frontal_cortex_ba9.csv","brain_hippocampus.csv","brain_hypothalamus.csv","brain_nucleus_accumbens_basal_ganglia.csv","brain_putamen_basal_ganglia.csv","brain_spinal_cord_cervical_c_1.csv","brain_substantia_nigra.csv","breast_mammary_tissue.csv","cells_ebv_transformed_lymphocytes.csv","cells_transformed_fibroblasts.csv","cervix_ectocervix.csv","cervix_endocervix.csv","colon_sigmoid.csv","colon_transverse.csv","esophagus_gastroesophageal_junction.csv","esophagus_mucosa.csv","esophagus_muscularis.csv","fallopian_tube.csv","heart_atrial_appendage.csv","heart_left_ventricle.csv","kidney_cortex.csv","liver.csv","lung.csv","minor_salivary_gland.csv","muscle_skeletal.csv","nerve_tibial.csv","pancreas.csv","pituitary.csv","prostate.csv","skin_not_sun_exposed_suprapubic.csv","skin_sun_exposed_lower_leg.csv","small_intestine_terminal_ileum.csv","spleen.csv","stomach.csv","thyroid.csv","uterus.csv","vagina.csv","whole_blood.csv")
feature_count = list()
testis_samples_expression <- read.table(files[1], header=TRUE, sep=",", na.strings=" ")[,1:6]
feature_count[1] = 5 # ncol(testis_samples_expression)-1
testis=testis_samples_expression[,1:ncol(testis_samples_expression)]
count.table <- data.frame(testis)
for (i in 2:length(files)){
  f = files[i]
  feature_count[f] = ncol(testis_samples_expression)-1
  other_tissue_expression <- read.table(f, header=TRUE, sep=",", na.strings=" ")[,1:6]
  feature_count[i] = 5 # ncol(testis_samples_expression)-1
  other=other_tissue_expression[,2:ncol(other_tissue_expression)]
  count.table <- cbind(count.table, other)
}
# other_tissue_expression <- read.table(files[2], header=TRUE, sep=",", na.strings=" ")
# feature_count[2] = ncol(other_tissue_expression)-1
# other=other_tissue_expression[,2:ncol(other_tissue_expression)]
# count.table <- cbind(count.table, other)
samplenames <- substring(files, 1, nchar(files)-4)
condition <- as.factor(rep(samplenames, feature_count))
col_names = colnames(count.table)[-1]
colData <- data.frame(condition = condition, row.names = col_names)
rownames(count.table)<-count.table[,1]
count.table <- count.table[-1]
dds<-DESeqDataSetFromMatrix(countData = count.table,colData = colData,design = ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds<-DESeq(dds)
# res_2 <- results(dds,contrast = c("condition",samplenames[1],samplenames[3]))
# diff_gene_deseq2 <-subset(res_2, padj < 0.05 & abs(log2FoldChange) > 4)
# out_result_file_name = paste(samplenames[1],samplenames[3], sep = "_")
# write.csv(diff_gene_deseq2,file= out_result_file_name)
setwd("F:/new_disease/3_label_properties/ovary_hyperexpression/results_fc2")
for (j in 2:length(files)){
  res <- results(dds,contrast = c("condition",samplenames[1],samplenames[j]))
  diff_gene_deseq2 <-subset(res, abs(log2FoldChange) > 1) #padj < 0.05 & 
  out_result_file_name = paste(samplenames[1],samplenames[j], sep = "_")
  write.csv(diff_gene_deseq2,file= out_result_file_name)
}

# Visualization (Volcano Map)

resultsNames(dds)

for (j in 2:length(files)){
  res <- results(dds,contrast = c("condition",samplenames[1],samplenames[j]))
  plotMA(res, main="DESeq2", ylim=c(-2,2))
}
tissue_name <- unlist(strsplit(samplenames[8], split = "_"))
tissue_name
res <- results(dds,contrast = c("condition",samplenames[1],samplenames[8]))
plotMA(res, main=paste(samplenames[1]," and ",tissue_name[1],tissue_name[2]), ylim=c(-2,2))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



