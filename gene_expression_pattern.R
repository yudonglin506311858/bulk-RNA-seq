setwd("D:/analysis/yangran")
data<-read.csv("terminalerythropoesis.csv",row.names = 1)
gene<-data
head(gene)
class(gene)

gene$mm_pro<-apply(gene[,1:4],1, mean, na.rm = T) 
gene$mm_baso<-apply(gene[,5:8],1, mean, na.rm = T) 
gene$mm_poly<-apply(gene[,9:13],1, mean, na.rm = T) 
gene$mm_ortho<-apply(gene[,14:17],1, mean, na.rm = T) 

gene <-gene[,17:20]
write.csv(gene,"四个时期的表达模式.CSV")

set.seed(1000)
gene <- scale(gene, center = TRUE, scale = TRUE)

library(ggplot2)

library(pheatmap)

library(reshape2)

cl <- kmeans(gene,9)

Cluster = paste0("cluster ",cl[["cluster"]]," : ",cl[["size"]][cl[["cluster"]]]," genes")

genenames <- rownames(gene)

gene1 <- data.frame(gene,Cluster,genenames)

gene2 <- melt(gene1)

# 绘制基因表达趋势折线图
pdf("表达聚类.pdf",height = 9,width = 18)
ggplot(gene2,aes(variable, value, group=genenames,col=Cluster)) + geom_line(size=2)+stat_summary(aes(group=1),fun=mean, geom="line", size=1.5, color="black") +
  
  facet_wrap(Cluster~.,scales = "free") +
  
  theme_bw() +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        
        axis.text = element_text(size=3, face = "bold"),
        
        strip.text = element_text(size = 15, face = "bold"),strip.background = element_rect(fill ="white", colour="white" ),panel.background = element_rect(fill="transparent"))+theme(panel.border = element_blank()) +theme(axis.line = element_line(size=1, colour = "black"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(panel.spacing.x = unit(1, "cm"),panel.spacing.y = unit(1, "cm"))+labs(y="expression",x="time")

dev.off()
#https://mp.weixin.qq.com/s/HnsZRjjbM6xSaAM1WJ8VpQ

pdf("表达模式.pdf",height = 9,width = 9)
table(gene2$Cluster)
#选择cluster1的基因
cluster1_gene<-gene2[gene2$Cluster=="cluster 1 : 509 genes","genenames"]
cluster1_gene<-as.data.frame(cluster1_gene)
cluster1_gene<-cluster1_gene[!duplicated(cluster1_gene),]
cluster1_gene<-as.data.frame(cluster1_gene)
write.csv(cluster1_gene,"cluster1_gene.csv")
#展示cluster2的基因的表达矩阵
cluster1_gene_expression<-gene[cluster1_gene$cluster1_gene,]
pheatmap::pheatmap(cluster1_gene_expression,cluster_rows = F,cluster_cols = F)


table(gene2$Cluster)
#选择cluster2的基因
cluster2_gene<-gene2[gene2$Cluster=="cluster 2 : 2622 genes","genenames"]
cluster2_gene<-as.data.frame(cluster2_gene)
cluster2_gene<-cluster2_gene[!duplicated(cluster2_gene),]
cluster2_gene<-as.data.frame(cluster2_gene)
write.csv(cluster2_gene,"cluster2_gene.csv")
#展示cluster2的基因的表达矩阵
cluster2_gene_expression<-gene[cluster2_gene$cluster2_gene,]
pheatmap::pheatmap(cluster2_gene_expression,cluster_rows = F,cluster_cols = F)

table(gene2$Cluster)

#选择cluster3的基因
cluster3_gene<-gene2[gene2$Cluster=="cluster 3 : 63 genes","genenames"]
cluster3_gene<-as.data.frame(cluster3_gene)
cluster3_gene<-cluster3_gene[!duplicated(cluster3_gene),]
cluster3_gene<-as.data.frame(cluster3_gene)
write.csv(cluster3_gene,"cluster3_gene.csv")
#展示cluster3的基因的表达矩阵
cluster3_gene_expression<-gene[cluster3_gene$cluster3_gene,]
pheatmap::pheatmap(cluster3_gene_expression,cluster_rows = F,cluster_cols = F)

table(gene2$Cluster)
#选择cluster4的基因
cluster4_gene<-gene2[gene2$Cluster=="cluster 4 : 31 genes","genenames"]
cluster4_gene<-as.data.frame(cluster4_gene)
cluster4_gene<-cluster4_gene[!duplicated(cluster4_gene),]
cluster4_gene<-as.data.frame(cluster4_gene)
write.csv(cluster4_gene,"cluster4_gene.csv")
#展示cluster4的基因的表达矩阵
cluster4_gene_expression<-gene[cluster4_gene$cluster4_gene,]
pheatmap::pheatmap(cluster4_gene_expression,cluster_rows = F,cluster_cols = F)


table(gene2$Cluster)
#选择cluster5的基因
cluster5_gene<-gene2[gene2$Cluster=="cluster 5 : 4 genes","genenames"]
cluster5_gene<-as.data.frame(cluster5_gene)
cluster5_gene<-cluster5_gene[!duplicated(cluster5_gene),]
cluster5_gene<-as.data.frame(cluster5_gene)
write.csv(cluster5_gene,"cluster5_gene.csv")
#展示cluster5的基因的表达矩阵
cluster5_gene_expression<-gene[cluster5_gene$cluster5_gene,]
pheatmap::pheatmap(cluster5_gene_expression,cluster_rows = F,cluster_cols = F)

table(gene2$Cluster)
#选择cluster6的基因
cluster6_gene<-gene2[gene2$Cluster=="cluster 6 : 8 genes","genenames"]
cluster6_gene<-as.data.frame(cluster6_gene)
cluster6_gene<-cluster6_gene[!duplicated(cluster6_gene),]
cluster6_gene<-as.data.frame(cluster6_gene)
write.csv(cluster6_gene,"cluster6_gene.csv")
#展示cluster6的基因的表达矩阵
cluster6_gene_expression<-gene[cluster6_gene$cluster6_gene,]
pheatmap::pheatmap(cluster6_gene_expression,cluster_rows = F,cluster_cols = F)


table(gene2$Cluster)
#选择cluster7的基因
cluster7_gene<-gene2[gene2$Cluster=="cluster 7 : 98 genes","genenames"]
cluster7_gene<-as.data.frame(cluster7_gene)
cluster7_gene<-cluster7_gene[!duplicated(cluster7_gene),]
cluster7_gene<-as.data.frame(cluster7_gene)
write.csv(cluster7_gene,"cluster7_gene.csv")
#展示cluster7的基因的表达矩阵
cluster7_gene_expression<-gene[cluster7_gene$cluster7_gene,]
pheatmap::pheatmap(cluster7_gene_expression,cluster_rows = F,cluster_cols = F)


table(gene2$Cluster)
#选择cluster8的基因
cluster8_gene<-gene2[gene2$Cluster=="cluster 8 : 5 genes","genenames"]
cluster8_gene<-as.data.frame(cluster8_gene)
cluster8_gene<-cluster8_gene[!duplicated(cluster8_gene),]
cluster8_gene<-as.data.frame(cluster8_gene)
write.csv(cluster8_gene,"cluster8_gene.csv")
#展示cluster8的基因的表达矩阵
cluster8_gene_expression<-gene[cluster8_gene$cluster8_gene,]
pheatmap::pheatmap(cluster8_gene_expression,cluster_rows = F,cluster_cols = F)

table(gene2$Cluster)
#选择cluster9的基因
cluster9_gene<-gene2[gene2$Cluster=="cluster 9 : 49547 genes","genenames"]
cluster9_gene<-as.data.frame(cluster9_gene)
cluster9_gene<-cluster9_gene[!duplicated(cluster9_gene),]
cluster9_gene<-as.data.frame(cluster9_gene)
write.csv(cluster9_gene,"cluster9_gene.csv")
#展示cluster9的基因的表达矩阵
cluster9_gene_expression<-gene[cluster9_gene$cluster9_gene,]
pheatmap::pheatmap(cluster9_gene_expression,cluster_rows = F,cluster_cols = F)

dev.off()





library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster1_gene$cluster1_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster1.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster1.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster1.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster1.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster1.csv")

pdf("KEGG_cluster1.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster1.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster2_gene$cluster2_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster2.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster2.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster2.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster2.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster2.csv")

pdf("KEGG_cluster2.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster2.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster3_gene$cluster3_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster3.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster3.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster3.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster3.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster3.csv")

pdf("KEGG_cluster3.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster3.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster4_gene$cluster4_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster4.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster4.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster4.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster4.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster4.csv")

pdf("KEGG_cluster4.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster4.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster5_gene$cluster5_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster5.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster5.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster5.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster5.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster5.csv")

pdf("KEGG_cluster5.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster5.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster6_gene$cluster6_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster6.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster6.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster6.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster6.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster6.csv")

pdf("KEGG_cluster6.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster6.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster7_gene$cluster7_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster7.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster7.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster7.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster7.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster7.csv")

pdf("KEGG_cluster7.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster7.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster8_gene$cluster8_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster8.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster8.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster8.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster8.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster8.csv")

pdf("KEGG_cluster8.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster8.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

b<-cluster9_gene$cluster9_gene
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_cluster9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_MF_cluster9.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_BP_cluster9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)

dev.off()

write.csv(ego,file="GO_BP_cluster9.csv")



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

barplot(ego)



dotplot(ego)



pdf("GO_CC_cluster9.pdf",width=10,height=6)

barplot(ego)
dotplot(ego)
dev.off()

write.csv(ego,file="GO_CC_cluster9.csv")



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_cluster9.csv")

pdf("KEGG_cluster9.pdf",width=10,height=6)

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_cluster9.pdf",width=10,height=6)
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



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


data<-read.csv("terminalerythropoesis.csv",row.names = 1)
normal<-data
colnames(normal)<-c("WT_pro","WT_pro","WT_pro","WT_pro","WT_baso","WT_baso","WT_baso","WT_baso","WT_poly","WT_poly","WT_poly","WT_poly","WT_ortho","WT_ortho","WT_ortho","WT_ortho")

normal$mm_pro<-apply(normal[,1:4],1, mean, na.rm = T) 
normal$mm_baso<-apply(normal[,5:8],1, mean, na.rm = T) 
normal$mm_poly<-apply(normal[,9:13],1, mean, na.rm = T) 
normal$mm_ortho<-apply(normal[,14:17],1, mean, na.rm = T) 

normal <-normal[,17:20]
normal <- data.matrix(normal)
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
pdf("mfuzz.pdf")
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER1.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER1.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER1.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER1.csv")

pdf("KEGG_CLUSTER1.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER1.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()




c2<-cl$cluster[cl$cluster == ]
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER2.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER2.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER2.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER2.csv")

pdf("KEGG_CLUSTER2.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER2.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



c3<-cl$cluster[cl$cluster == ]
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER3.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER3.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER3.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER3.csv")

pdf("KEGG_CLUSTER3.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER3.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER4.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER4.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER4.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER4.csv")

pdf("KEGG_CLUSTER4.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER4.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER5.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER5.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER5.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER5.csv")

pdf("KEGG_CLUSTER5.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER5.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER6.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER6.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER6.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER6.csv")

pdf("KEGG_CLUSTER6.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER6.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER7.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER7.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER7.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER7.csv")

pdf("KEGG_CLUSTER7.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER7.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER8.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER8.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER8.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER8.csv")

pdf("KEGG_CLUSTER8.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER8.pdf")
barplot(ego)
dotplot(ego)
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

barplot(ego)



dotplot(ego)



pdf("GO_MF_CLUSTER9.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_BP_CLUSTER9.pdf")

barplot(ego)
dotplot(ego)

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

barplot(ego)



dotplot(ego)



pdf("GO_CC_CLUSTER9.pdf")

barplot(ego)
dotplot(ego)
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

barplot(ekegg)



dotplot(ekegg)

write.csv(ekegg,file="kegg_CLUSTER9.csv")

pdf("KEGG_CLUSTER9.pdf")

barplot(ekegg)

dotplot(ekegg)


dev.off()

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_CLUSTER9.pdf")
barplot(ego)
dotplot(ego)
#柱状图,各自展示前8个term
barplot(ego,drop = TRUE, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =8,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



