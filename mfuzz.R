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
set.seed(123)
library(clusterProfiler)
library(org.Mm.eg.db)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
## 加载R包
library("Mfuzz")
rm(list=ls())
getwd()
setwd("D:/anxiiuli_data/count_炎红/CD44")


cts<-read.csv("ternimalE_fpkm.csv",row.names = 1)

#读取tpm文件：
head(cts)


data<-cts[,c(6:16)]
normal<-data
head(normal)

normal$day0<-apply(normal[,1:3],1, mean, na.rm = T) 
normal$day4<-apply(normal[,4:5],1, mean, na.rm = T) 
normal$day6<-apply(normal[,6:8],1, mean, na.rm = T) 
normal$day8<-apply(normal[,9:11],1, mean, na.rm = T) 

normal <-normal[,12:15]
normal <- data.matrix(normal)
write.csv(normal,"mfuzz_表达矩阵.csv")
eset <- new("ExpressionSet",exprs = normal)
## 过滤缺失超过25%的基因
gene.r <- filter.NA(eset, thres=0.25)
## mean填补缺失
gene.f <- fill.NA(gene.r,mode="mean")
## knn/wknn方法表现更好，但是计算起来比较复杂
gene.f <- fill.NA(gene.r,mode="knn")
#gene.f <- fill.NA(gene.r,mode="wknn")
## 过滤标准差为0的基因
tmp <- filter.std(gene.f,min.std=0)
## 标准化
gene.s <- standardise(tmp)
## 聚类个数
c <- 9
## 计算最佳的m值
m <- mestimate(gene.s)
## 聚类
cl <- mfuzz(gene.s, c = c, m = m)
## 查看每类基因数目
cl$size
## 查看每类基因ID
cl$cluster[cl$cluster == 1]
## 输出基因ID
write.table(cl$cluster,"output.txt",quote=F,row.names=T,col.names=F,sep="\t")
## 绘制折线图
pdf("mfuzz.pdf",width=20,height=20)
mfuzz.plot(gene.s,cl,mfrow=c(3,3),new.window= FALSE)
dev.off()
c1<-cl$cluster[cl$cluster == 1]
c1 <- as.data.frame(c1)
c1<-rownames(c1)
write.csv(c1,"cluster1_gene.csv")
c2<-cl$cluster[cl$cluster == 2]
c2 <- as.data.frame(c2)
c2<-rownames(c2)
write.csv(c2,"cluster2_gene.csv")
c3<-cl$cluster[cl$cluster == 3]
c3 <- as.data.frame(c3)
c3<-rownames(c3)
write.csv(c3,"cluster3_gene.csv")
c4<-cl$cluster[cl$cluster == 4]
c4 <- as.data.frame(c4)
c4<-rownames(c4)
write.csv(c4,"cluster4_gene.csv")
c5<-cl$cluster[cl$cluster == 5]
c5 <- as.data.frame(c5)
c5<-rownames(c5)
write.csv(c5,"cluster5_gene.csv")
c6<-cl$cluster[cl$cluster == 6]
c6 <- as.data.frame(c6)
c6<-rownames(c6)
write.csv(c6,"cluster6_gene.csv")
c7<-cl$cluster[cl$cluster == 7]
c7 <- as.data.frame(c7)
c7<-rownames(c7)
write.csv(c7,"cluster7_gene.csv")
c8<-cl$cluster[cl$cluster == 8]
c8 <- as.data.frame(c8)
c8<-rownames(c8)
write.csv(c8,"cluster8_gene.csv")
c9<-cl$cluster[cl$cluster == 9]
c9 <- as.data.frame(c9)
c9<-rownames(c9)
write.csv(c9,"cluster9_gene.csv")


dir.create("mfuzz_Enrichment")
setwd("./mfuzz_Enrichment")


set.seed(123)
library(clusterProfiler)
library(org.Mm.eg.db)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
c1<-cl$cluster[cl$cluster == 1]
c1 <- as.data.frame(c1)
c1<-rownames(c1)
b<-c1
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER1.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER1.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER1.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER1.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER1.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER1.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER1.csv")

pdf("KEGG_CLUSTER1.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER1.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




c2<-cl$cluster[cl$cluster == 2]
c2 <- as.data.frame(c2)
c2<-rownames(c2)
b<-c2
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER2.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER2.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER2.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER2.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER2.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER2.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER2.csv")

pdf("KEGG_CLUSTER2.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER2.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



c3<-cl$cluster[cl$cluster == 3]
c3 <- as.data.frame(c3)
c3<-rownames(c3)
b<-c3
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER3.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER3.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER3.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER3.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER3.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER3.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER3.csv")

pdf("KEGG_CLUSTER3.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER3.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




c4<-cl$cluster[cl$cluster ==4 ]
c4 <- as.data.frame(c4)
c4<-rownames(c4)
b<-c4
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER4.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER4.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER4.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER4.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER4.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER4.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER4.csv")

pdf("KEGG_CLUSTER4.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER4.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




c5<-cl$cluster[cl$cluster == 5]
c5 <- as.data.frame(c5)
c5<-rownames(c5)
b<-c5
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER5.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER5.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER5.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER5.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER5.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER5.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER5.csv")

pdf("KEGG_CLUSTER5.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER5.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()





c6<-cl$cluster[cl$cluster == 6]
c6 <- as.data.frame(c6)
c6<-rownames(c6)
b<-c6
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER6.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER6.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER6.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER6.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER6.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER6.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER6.csv")

pdf("KEGG_CLUSTER6.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER6.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


c7<-cl$cluster[cl$cluster == 7]
c7 <- as.data.frame(c7)
c7<-rownames(c7)
b<-c7
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER7.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER7.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER7.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER7.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER7.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER7.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER7.csv")

pdf("KEGG_CLUSTER7.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER7.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

c8<-cl$cluster[cl$cluster == 8]
c8 <- as.data.frame(c8)
c8<-rownames(c8)
b<-c8
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER8.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER8.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER8.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER8.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER8.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER8.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER8.csv")

pdf("KEGG_CLUSTER8.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER8.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



c9<-cl$cluster[cl$cluster == 9]
c9 <- as.data.frame(c9)
c9<-rownames(c9)
b<-c9
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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



dotplot(ego, showCategory =50)



pdf("GO_MF_CLUSTER9.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_MF_CLUSTER9.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_BP_CLUSTER9.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)

dev.off()

write.csv(ego,file="GO_BP_CLUSTER9.csv")



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



dotplot(ego, showCategory =50)



pdf("GO_CC_CLUSTER9.pdf",width=20,height=20)

barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()

write.csv(ego,file="GO_CC_CLUSTER9.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg, showCategory =50)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER9.csv")

pdf("KEGG_CLUSTER9.pdf",width=20,height=20)

barplot(ekegg, showCategory =50)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER9.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



