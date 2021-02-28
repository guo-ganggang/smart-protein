install.packages('stringi')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('farver')
install.packages('labeling')


windowsFonts(Times=windowsFont("Times New Roman"))
library(stringi)
library(ggplot2)
library(dplyr)
library(farver)
library(labeling)

##Read in and collate KEGG results

downgokegg<-read.delim('KEGG_PATHWAY.txt')

enrich<-downgokegg

enrich_signif=enrich[which(enrich$PValue<0.05),]

enrich_signif=enrich_signif[,c(1:3,5)]

head(enrich_signif)
enrich_signif=data.frame(enrich_signif)

KEGG=enrich_signif

KEGG$Term<-stri_sub(KEGG$Term,10,100)

kegg_ea <- ggplot(KEGG,aes(x=Count,y=Term))+geom_point(aes(color=PValue,size=Count))+scale_color_gradient(low='green',high='red')+theme_bw()+theme(text = element_text(size=12),axis.text.x = element_text(family="Times",face = "bold",color = "black",size=10),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(kegg_ea, file='kegg_ea.pdf', width=12, height=6)


###The results of GO enrichment were read and sorted

GO_CC<-read.delim('GOTERM_CC_DIRECT.txt')

GO_CC_signif=GO_CC[which(GO_CC$PValue<0.05),]

GO_CC_signif=GO_CC[,c(1:3,5)]

head(GO_CC_signif)

GO_CC_signif=data.frame(GO_CC_signif)

GO_CC_signif$Term<-stri_sub(GO_CC_signif$Term,12,100)

GO_BP<-read.delim('GOTERM_BP_DIRECT.txt')

GO_BP_signif=GO_BP[which(GO_BP$PValue<0.05),]

GO_BP_signif=GO_BP_signif[,c(1:3,5)]

head(GO_BP_signif)

GO_BP_signif=data.frame(GO_BP_signif)

GO_BP_signif$Term<-stri_sub(GO_BP_signif$Term,12,100)

GO_MF<-read.delim('GOTERM_MF_DIRECT.txt')

GO_MF_signif=GO_MF[which(GO_MF$PValue<0.05),]

GO_MF_signif=GO_MF_signif[,c(1:3,5)]

head(GO_MF_signif)

GO_MF_signif=data.frame(GO_MF_signif)

GO_MF_signif$Term<-stri_sub(GO_MF_signif$Term,12,100)

enrich_signif=rbind(GO_BP_signif,rbind(GO_CC_signif,GO_MF_signif))

go=enrich_signif

go=arrange(go,go$Category,go$PValue)

##Legend name setting

m=go$Category

m=gsub("TERM","",m)

m=gsub("_DIRECT","",m)

go$Category=m

GO_term_order=factor(as.integer(rownames(go)),labels = go$Term)

COLS<-c("#66C3A5","#8DA1CB","#FD8D62")

###Began to draw
###labels = c("Biological process", "Cellular component", "Molecular function")

p <- ggplot(data=go,aes(x=GO_term_order,y=Count,fill=Category))+geom_bar(stat = "identity",width = 0.8)+scale_fill_manual(values = COLS,labels = c("Biological process", "Cellular component", "Molecular function"))+theme_bw()+theme(text=element_text(family="Times", face="bold", size=12))+xlab("GO terms")+ylab("Gene counts")+labs()+theme(text = element_text(size=12),axis.text.x = element_text(family="Times",face = "bold",color = "black",angle = 90,vjust = 1,hjust = 1,size=12))
go_ea <- p + theme(legend.background = element_rect(fill="aliceblue", size=0.5, linetype="solid")) + theme(legend.position = c(0.9, 0.8),legend.text=element_text(family="Times",color = "black",size=12))
ggsave(go_ea, file='go_ea.pdf', width=12, height=10)

