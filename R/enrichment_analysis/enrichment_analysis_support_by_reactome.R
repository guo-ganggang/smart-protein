
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")


library(ReactomePA)
library(ggplot2)


##Loading the data
df <- read.table("predicted_causal_gene_ids.txt",header=T,quote = "",sep = "\t",colClasses = c("character","character","character","character"))
de <- df[c("To")]
ids=na.omit(de$To)  #提取出非NA的ENTREZID
genes=ids
head(genes)

##enrichment analysis
epw <-enrichPathway(gene=genes,pvalueCutoff=0.05, readable=T)
head(as.data.frame(epw))

##bubble chart
write.table(as.data.frame(epw@result), file="enrichPathway_results.txt")
p = dotplot(epw,showCategory = 15,font.size=12)+ ggplot2::xlab("Gene ratio")+
  theme_bw() +
  theme(text=element_text(family="Times", face="bold", size=12))+theme(text = element_text(size=12),axis.text.x = element_text(family="Times",face = "bold",color = "black",size=12))
ggsave(p, filename = "kegg_ea.pdf", device = "pdf")


