Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(pheatmap)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
setwd("D:/anxiiuli_data/count_炎红/")

setwd("D:/anxiiuli_data/count_炎红/CD44")




#读入数据并合并
non_immune_ery13<-read.table("13-非炎红_FRAS220042674-1r.count",header=T)
non_immune_ery14<-read.table("14-非炎红_FRAS220042675-1r.count",header=T)
non_immune_ery15<-read.table("15-非炎红_FRAS220042676-1r.count",header=T)
#0天造模，野生型
immune_ery13<-read.table("13-炎红_FRAS220042667-1r.count",header=T)
immune_ery14<-read.table("14-炎红_FRAS220042668-1r.count",header=T)
immune_ery15<-read.table("15-炎红_FRAS220042669-1r.count",header=T)
#4天造模
wt_phz_3<-read.table("3-炎红_FRAS220042670-1r.count",header=T)
wt_phz_5<-read.table("5-炎红_FRAS220042671-1r.count",header=T)
#6天造模
wt_phz_7<-read.table("7-炎红_FRAS220042661-2r.count",header=T)
wt_phz_8<-read.table("8-炎红_FRAS220042662-2r.count",header=T)
wt_phz_9<-read.table("9-炎红_FRAS220042663-2r.count",header=T)
#8天造模
wt_phz_10<-read.table("10-炎红_FRAS220042664-2r.count",header=T)
wt_phz_11<-read.table("11-炎红_FRAS220042665-2r.count",header=T)
wt_phz_12<-read.table("12-炎红_FRAS220042666-1r.count",header=T)
#CD44 KO 4天造模
ko_phz_9<-read.table("9-炎红黑_FRAS220042672-1r.count",header=T)
ko_phz_10<-read.table("10-炎红黑_FRAS220042673-1r.count",header=T)



#合并数据
data<-merge(non_immune_ery13,non_immune_ery14,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,non_immune_ery15,by=c("Geneid","Chr","Start","End","Strand","Length"))
#0天造模，野生型
data<-merge(data,immune_ery13,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,immune_ery14,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,immune_ery15,by=c("Geneid","Chr","Start","End","Strand","Length"))
#4天造模
data<-merge(data,wt_phz_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,wt_phz_5,by=c("Geneid","Chr","Start","End","Strand","Length"))
#6天造模
data<-merge(data,wt_phz_7,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,wt_phz_8,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,wt_phz_9,by=c("Geneid","Chr","Start","End","Strand","Length"))
#8天造模
data<-merge(data,wt_phz_10,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,wt_phz_11,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,wt_phz_12,by=c("Geneid","Chr","Start","End","Strand","Length"))
#CD44 KO 4天造模
data<-merge(data,ko_phz_9,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,ko_phz_10,by=c("Geneid","Chr","Start","End","Strand","Length"))



WT_count<-data

head(WT_count)
dim(WT_count)
colnames(WT_count)


#删除ensmusg的版本号

#WT_count$Geneid<-gsub("\\.*","",WT_count$Geneid)
#WT_count$Geneid <- gsub("\\.[0-9]*$", "", WT_count$Geneid)
rownames(WT_count)<-WT_count$Geneid
#https://www.nhooo.com/note/qa02b6.html
#https://www.biostars.org/p/178726/
WT_count_all<-WT_count[,c(1,6:22)]
head(WT_count_all)
colnames(WT_count_all)<-c("Geneid","Length","non_immune_ery_1","non_immune_ery_2","non_immune_ery_3",
                          "immune_ery_day0_1","immune_ery_day0_2","immune_ery_day0_3",
                          "wt_phz_immune_ery_day4_1","wt_phz_immune_ery_day4_2",
                          "wt_phz_immune_ery_day6_1","wt_phz_immune_ery_day6_2","wt_phz_immune_ery_day6_3",
                          "wt_phz_immune_ery_day8_1","wt_phz_immune_ery_day8_2","wt_phz_immune_ery_day8_3",
                          "ko_phz_immune_ery_day4_1","ko_phz_immune_ery_day4_2"
)


cts<-WT_count_all
head(cts)
dim(cts)

write.csv(cts,"ensembl_gene_expression_count_yan.csv")



#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("mmapiens_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_human.rds")
mart<-readRDS("mart_human.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"human_gene_ensembl_transition.csv")




library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
write.csv(cts,file = "ternimalE_count_genesymbol.csv",row.names = F)
cts<-cts[,-1]
cts$external_gene_name<-make.names(cts$external_gene_name, unique = TRUE)
rownames(cts)<-cts$external_gene_name
write.table(cts,"data.txt", sep = "\t")



#读取count文件：
write.csv(cts,file="ternimalE_count.csv")

#读取tpm文件：
head(cts)
kb <- cts$Length / 1000
head(kb)
countdata <- cts[,3:18]
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm,file="ternimalE_tpm.csv")

#读取fpkm文件：
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm,file="ternimalE_fpkm.csv")


#PCA
condition<-factor(c("non_immune_ery","non_immune_ery","non_immune_ery",
                    "immune_ery_day0","immune_ery_day0","immune_ery_day0",
                    "wt_phz_immune_ery_day4","wt_phz_immune_ery_day4",
                    "wt_phz_immune_ery_day6","wt_phz_immune_ery_day6","wt_phz_immune_ery_day6",
                    "wt_phz_immune_ery_day8","wt_phz_immune_ery_day8","wt_phz_immune_ery_day8",
                    "ko_phz_immune_ery_day4","ko_phz_immune_ery_day4"
),
                  levels = c("non_immune_ery","immune_ery_day0","wt_phz_immune_ery_day4","wt_phz_immune_ery_day6","wt_phz_immune_ery_day8","ko_phz_immune_ery_day4"))
head(cts)
tmp<-cts[,3:18]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:18]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:18],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")
vsd_all$condition<-  factor(c("non_immune_ery_1","non_immune_ery_2","non_immune_ery_3",
    "immune_ery_day0_1","immune_ery_day0_2","immune_ery_day0_3",
    "wt_phz_immune_ery_day4_1","wt_phz_immune_ery_day4_2",
    "wt_phz_immune_ery_day6_1","wt_phz_immune_ery_day6_2","wt_phz_immune_ery_day6_3",
    "wt_phz_immune_ery_day8_1","wt_phz_immune_ery_day8_2","wt_phz_immune_ery_day8_3",
    "ko_phz_immune_ery_day4_1","ko_phz_immune_ery_day4_2"
  ),levels = c("non_immune_ery_1","non_immune_ery_2","non_immune_ery_3",
               "immune_ery_day0_1","immune_ery_day0_2","immune_ery_day0_3",
               "wt_phz_immune_ery_day4_1","wt_phz_immune_ery_day4_2",
               "wt_phz_immune_ery_day6_1","wt_phz_immune_ery_day6_2","wt_phz_immune_ery_day6_3",
               "wt_phz_immune_ery_day8_1","wt_phz_immune_ery_day8_2","wt_phz_immune_ery_day8_3",
               "ko_phz_immune_ery_day4_1","ko_phz_immune_ery_day4_2"
  ))
# 
#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

cts<-cts[,3:6]
cor(cts,method = "pearson")
# ctrl_1    ctrl_2  ptp4a3_1  ptp4a3_2
# ctrl_1   1.0000000 0.9863111 0.8094934 0.5713300
# ctrl_2   0.9863111 1.0000000 0.7706288 0.5144732
# ptp4a3_1 0.8094934 0.7706288 1.0000000 0.9099334
# ptp4a3_2 0.5713300 0.5144732 0.9099334 1.0000000
pheatmap(cor(cts))

pheatmap(cor(cts[,3:18],method = "pearson"),display_numbers = T)
#pheatmap(cor(cts,method = "spearman"),display_numbers = T)

#做一下差异表达分析
res_all <- results(dds_all)
head(res_all)
res_all <- res_all[order(res_all$padj),]
head(res_all)
diff_gene <- subset(res_all, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res_all)
write.csv(diff_gene,file = "diff_gene_all.csv",row.names = T)
diff_gene_up <- subset(res_all, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up.csv",row.names = T)
diff_gene_down <- subset(res_all, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down.csv",row.names = T)
head(diff_gene_down)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)

resdata <- merge(as.data.frame(res_all), as.data.frame(counts(dds_all, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all.csv",row.names = F)
head(resdata)
resdata = res_all[order(resdata$pvalue),]
summary(resdata)
table(resdata$pvalue<0.05)#number of true 小于0.05 的基因个数

diff_gene<-as.data.frame(diff_gene)
term<-rownames(diff_gene)
newdata<-cts[,3:18][c(term),]

pheatmap(newdata,scale ="row",border_color = NA)
pheatmap(newdata,scale ="row",border_color = NA,show_rownames = F)
pheatmap(newdata,scale ="row",border_color = NA,filename = "heatmap.pdf",width = 7,height = 300)




# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread("DEG_all.csv")

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
labelGene <- rbind(labelGene1,labelGene2)
#labelGene <- rbind(labelGene,diffData[1685,])#Ptp4a3

options(ggrepel.max.overlaps = Inf)
#pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())+labs(title= "vocano plot")

ggsave("volcano.pdf",width = 6,height = 6)


cts<-read.csv("ternimalE_count.csv",row.names = 1)

#读取tpm文件：
head(cts)
data<-cts
data<-data[,-1]
data[1:4,1:4]#就是普通的表达矩阵
dim(data)
data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

class(data)
group=c("non_immune_ery","immune_ery_day0","wt_phz_immune_ery_day4","wt_phz_immune_ery_day6","wt_phz_immune_ery_day8","ko_phz_immune_ery_day4")#就是分成多少组
num=c(3,3,2,3,3,2) #这几组每组分别有多少个
num
dir.create("Batch_Enrichment")
setwd("./Batch_Enrichment")



Batch_Enrichment<-function(exprSet,group,num){
  ## creat a group
  group_list= factor(rep(group,num))
  group_list
  #exprSet<-data
  colData=data.frame(row.names = colnames(exprSet),
                     group=group_list)
  
  
  for (i in 1:length(group)){
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                                       colData = colData,
                                       design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
      resname=DESeq2::resultsNames(dds)[j]
      
      res=DESeq2::results(dds, name=resname)
      print(resname)
      res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      res_lfc
      res = res[order(res$pvalue),]
      res<-na.omit(as.data.frame(res))
      print(head(res))
      summary(res)
      
      res_lfc = res_lfc[order(res_lfc$pvalue),]
      res_lfc<-na.omit(as.data.frame(res_lfc))
      print(head(res_lfc))
      summary(res_lfc)
      
      res<-res_lfc
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      #write.csv(resdata,paste0(resname,".csv"),row.names = F)
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      resdata <- merge(as.data.frame(res_lfc), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(resdata,paste0(resname,".csv"),row.names = F)
      
      head(resdata)
      summary(res[order(resdata$pvalue),])
      table(resdata$padj<0.05)#number of true 小于0.05 的基因个数
      diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
      head(diff_gene)
      dim(diff_gene)
      diff_gene <- merge(as.data.frame(diff_gene), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene,file = paste0("diff_gene_all_",resname,".csv"),row.names = F)
      diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
      diff_gene_up <- merge(as.data.frame(diff_gene_up), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_up,file = paste0("diff_gene_up_",resname,".csv"),row.names = F)
      diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
      diff_gene_down <- merge(as.data.frame(diff_gene_down), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_down,file = paste0("diff_gene_down_",resname,".csv"),row.names = F)
      
      #pdf(paste0(resname,".pdf"))
      library(ggplot2)
      dataset <-read.csv(paste0(resname,".csv"),header = TRUE)
      dim(dataset)
      head(dataset)
      dataset <-na.omit(dataset)
      cut_off_pvalue = 0.01
      cut_off_log2FoldChange = 1
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= paste0(resname,"_volcano plot")) +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )
      
      
      #options(ggrepel.max.overlaps = Inf)
      library(ggrepel)
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      
      dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -1| dataset$log2FoldChange >= 1,as.character(dataset$Row.names),"")
      
      p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)
      ggsave(paste0(resname,"_volcano1.pdf"))
      
      
      # grammar
      library(tidyverse)
      library(magrittr)
      library(glue)
      library(data.table)
      
      # analysis
      library(DESeq2)
      
      # graphics
      library(ggplot2)
      library(ggrepel)
      library(ggsci)
      library(scales)
      library(latex2exp)
      diffData <- fread(paste0(resname,".csv"))
      #diffData <- fread("DEG_all.csv")
      colnames(diffData)[1] <- "gene"
      
      diffData[is.na(padj), padj := 1][]
      diffData[, p := -log10(padj)][]
      
      
      diffData[, type := "ns"][]
      diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]
      
      labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
      labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
      labelGene <- rbind(labelGene1,labelGene2)
      
      #pal_nejm()(8) %>% show_col()
      typeColor <- structure(
        c(pal_nejm()(2), "gray80"),
        names = c("up", "down", "ns")
      )
      
      ggplot(diffData, aes(x = log2FoldChange, y = p)) +
        geom_point(aes(color = type, size = p), show.legend = F) +
        geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
        geom_text_repel(
          data = labelGene, aes(label = gene),
          size = 3, fontface = 3,
          nudge_x = .5, nudge_y = .5) +
        scale_radius(range = c(.1, 2)) +
        scale_color_manual(values = typeColor) +
        scale_y_continuous(expand = expansion(c(0, 0.05))) +
        labs(
          x = TeX("$log_{2}(Fold\\,Change)$"),
          y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
        theme(
          aspect.ratio = 1,
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line())+labs(title= paste0(resname,"_vocano plot"))
      ggsave(paste0(resname,"_volcano2.pdf"))
      
      try({
   
        
        #下调基因
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        diff_gene_down<-read.csv(file = paste0("diff_gene_down_",resname,".csv"))
        
        b<-as.data.frame(diff_gene_down$Row.names)
        eg = bitr(b$`diff_gene_down$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_down_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_down_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_down_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        library(DOSE)
        ekegg1<-setReadable(ekegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
        barplot(ekegg1, showCategory =50)
        ggsave(paste0("barplot_KEGG_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg1, showCategory =50)
        ggsave(paste0("dotplot_KEGG_down_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg1,file=paste0("KEGG_down_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_down_",resname,".csv"))
        # 
      })
        
        try({
        #上调基因
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        diff_gene_up<-read.csv(file = paste0("diff_gene_up_",resname,".csv"))
        
        b<-as.data.frame(diff_gene_up$Row.names)
        eg = bitr(b$`diff_gene_up$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_up_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_up_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_up_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        
        library(DOSE)
        ekegg1<-setReadable(ekegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
        
        barplot(ekegg1, showCategory =50)
        ggsave(paste0("barplot_KEGG_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg1, showCategory =50)
        ggsave(paste0("dotplot_KEGG_up_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg1,file=paste0("KEGG_up_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_up_",resname,".csv"))
        })
        try({
        # 
        # #不区分上下调基因：all
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        library(DO.db)
        diff_gene<-read.csv(file = paste0("diff_gene_all_",resname,".csv"))
        
        b<-as.data.frame(diff_gene$Row.names)
        eg = bitr(b$`diff_gene$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_all_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_all_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_all_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        
        library(DOSE)
        ekegg1<-setReadable(ekegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
        
        barplot(ekegg1, showCategory =50)
        ggsave(paste0("barplot_KEGG_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg1, showCategory =50)
        ggsave(paste0("dotplot_KEGG_all_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg1,file=paste0("KEGG_all_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_all_",resname,".csv"))
        # 
        
      })  
    }
    
  }
}
Batch_Enrichment(data,group,num)


