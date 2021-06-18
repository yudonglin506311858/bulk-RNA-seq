rm(list = ls())
options(stringsAsFactors = F)
setwd("H:/Codes_for_Fang_Dong_et_al_NCB_2020-master/new")

# 首先读取两个文件
sig_matrix <-"H:/Codes_for_Fang_Dong_et_al_NCB_2020-master/new/LM22.txt"  # CIBERSORT 内置数据库挖掘
mixture_file <- "H:/Codes_for_Fang_Dong_et_al_NCB_2020-master/mixture_file.txt"  # 约80M，TCGA数据库

# 两个表达矩阵需要取交集
LM22.file <- sig_matrix
exp.file <- mixture_file
#read in data
X <- read.csv(LM22.file,row.names=1,sep="\t")
Y <- read.table(exp.file, header=T, sep="\t", check.names=F)
#Y <- Y[!duplicated(Y[,1]),]###去重复基因名
#rownames(Y)<-Y[,1]
#Y<-Y[,-1]
X <- data.matrix(X)###convert data as matrix
Y <- data.matrix(Y)
Y[1:4,1:4]##check data
X[1:4,1:4]
dim(X)
dim(Y)

X <- X[order(rownames(X)),]###行名字母排序
Y <- Y[order(rownames(Y)),]###行名字母排序

#anti-log if max < 50 in mixture file
if(max(Y) < 50) {Y <- 2^Y} ###如果Y矩阵中最大值<50，则变为2的y次方，也就是原始Y是被log2的

QN = F #QN = Quantile normalization of input mixture (default = TRUE)
#quantile normalization of mixture file
if(QN == TRUE){
  tmpc <- colnames(Y)
  tmpr <- rownames(Y)
  Y <- normalize.quantiles(Y)#preprocessCore的函数，正态化数据
  colnames(Y) <- tmpc
  rownames(Y) <- tmpr
}

#intersect genes
Xgns <- row.names(X)
Ygns <- row.names(Y)
YintX <- Ygns %in% Xgns ###y中取x
Y <- Y[YintX,] ###取共有子集
XintY <- Xgns %in% row.names(Y)###x中取y
X <- X[XintY,]
dim(X)
dim(Y)


#standardize sig matrix
X <- (X - mean(X)) / sd(as.vector(X)) ###标准化数据
Y[1:4,1:4]
X[1:4,1:4]
boxplot(X[,1:4])
save(X,Y,file = 'input.Rdata')




rm(list = ls())
options(stringsAsFactors = F)
load(file = 'input.Rdata')
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)
# 下面的演示是为了搞清楚 CoreAlg 函数
# 并不需要保存任何信息

# 从表达矩阵Y里面，随机挑选LM22矩阵基因数量的表达量值
Ylist <- as.list(data.matrix(Y)) ###将Y矩阵中每一个数值作为list的一个元素
yr <- as.numeric(Ylist[ sample(length(Ylist),dim(X)[1]) ])###从Ylist随机挑选nrow（x）个元素
# yr 这个时候是一个假设的样本
#standardize mixture,就是scale 函数
yr <- (yr - mean(yr)) / sd(yr) ##标准化样本
boxplot(yr)
# 每次随机挑选的yr，都是需要走后面的流程

# 一切都是默认值的支持向量机
# 这里的X是LM22矩阵，不同的免疫细胞比例组合成为不同的yr
# 这里的yr是随机的，反推免疫细胞比例
library(e1071)
out=svm(X,yr)##e1071包的函数
out
out$SV ##SV=支持向量

# 需要修改的参数包括：type="nu-regression",kernel="linear",nu=nus,scale=F
###参数含义，采用nu-regression,nu回归，采用线性方法从训练数据中学习得到模型，nu参数为nus值，不scale
###SVM参数https://stats.stackexchange.com/questions/94118/difference-between-ep-svr-and-nu-svr-and-least-squares-svr
###https://blog.csdn.net/csqazwsxedc/article/details/52230092

svn_itor <- 3

y=yr
#try different values of nu
res <- function(i){
  if(i==1){nus <- 0.25}
  if(i==2){nus <- 0.5}
  if(i==3){nus <- 0.75}
  model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
  model
}

#Execute In a parallel way the SVM
####Windows没有办法用mclapply开多核的，可以用parlapply
library(parallel)

if(Sys.info()['sysname'] == 'Windows') {
  
  out <- mclapply(1:svn_itor, res, mc.cores=1) 
}else {
  out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
}
# 运行了Support Vector Machines，函数是 svm {e1071}

###windows开多核
library(parallel)
clnum<-detectCores() 
cl <- makeCluster(getOption("cl.cores", clnum))
clusterExport(cl, c("X","y"),envir=environment())
clusterEvalQ(cl,library(e1071))
out <- parLapply(cl,1:svn_itor,res)
####
out
#Initiate two variables with 0
nusvm <- rep(0,svn_itor)
corrv <- rep(0,svn_itor)


t <- 1
while(t <= svn_itor) {
  
  # 得到两个向量之间矩阵乘法的权重，此时应该只得到一个数字。
  # 这样做是乘以系数
  
  # 支持向量是数据集的点，它们靠近分隔类别的平面
  # 现在的问题是，我没有任何类别（离散变量，例如“运动”、“电影”），但我有一个连续变量
  mySupportVectors <- out[[t]]$SV ###不同nu参数的支持向量
  
  # 系数定义
  myCoefficients <- out[[t]]$coefs ###不同nu参数的系数
  weights = t(myCoefficients) %*% mySupportVectors ####2个矩阵的乘积
  
  # 设置权重和相关性
  weights[which(weights<0)]<-0 ##小于0的乘积为0
  w<-weights/sum(weights) ##相关性
  
  # 根据对应的权重与参考集相乘
  u <- sweep(X,MARGIN=2,w,'*') ###sweep类似于apply，多了一个STATS，代表是运算的参数
  
  # 统计每行总和
  k <- apply(u, 1, sum)
  nusvm[t] <- sqrt((mean((k - y)^2))) ###标准差
  corrv[t] <- cor(k, y) ##相关性
  t <- t + 1 ###t从1开始循环，直到t=3
}
#pick best model
rmses <- nusvm
corrv
mn <- which.min(rmses) ###去标准差最小的nu值为best model
mn  
#[1] 1
model <- out[[mn]]
# 从nus为0.25,0.5,0.75的3个模型里面挑选一个即可

#get and normalize coefficients

q <- t(model$coefs) %*% model$SV 

q[which(q<0)]<-0
# w 就是计算后的22种免疫细胞的比例
w <- (q/sum(q))

mix_rmse <- rmses[mn]
mix_r <- corrv[mn]

# 会返回这个随机的y的免疫细胞组成情况，就是权重w
newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
newList

# 根据对应的权重与参考集相乘
u <- sweep(X,MARGIN=2,w,'*') 
k <- apply(u, 1, sum)
plot(y,k)
sqrt((mean((k - y)^2))) 
cor(k, y)
# 通常这个预测结果都惨不忍睹





# 每次把表达矩阵通过去卷积拆分成为LM22的免疫细胞比例结果

# 并且包装成为函数，如下：

#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#Core algorithm
CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

library(preprocessCore)
library(parallel)
library(e1071)
source("cibersort.R") 

itor <- 1
Ylist <- as.list(data.matrix(Y))
dist <- matrix()
# 就是把 CoreAlg 函数运行1000次
perm=1000
while(itor <= perm){
  print(itor) # 打印进度
  
  #random mixture
  yr <- as.numeric(Ylist[ sample(length(Ylist),dim(X)[1]) ])
  
  #standardize mixture
  yr <- (yr - mean(yr)) / sd(yr)
  
  #run CIBERSORT core algorithm
  result <- CoreAlg(X, yr)
  
  mix_r <- result$mix_r
  
  #store correlation
  if(itor == 1) {dist <- mix_r}
  else {dist <- rbind(dist, mix_r)}
  
  itor <- itor + 1
}
####
newList <- list("dist" = dist)###获取1000次计算相关系数
nulldist=sort(newList$dist) ###w值排序
# 这个nulldist 主要是用来计算P值
if(F){
  
  P=perm
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)} ###取最小的非负p值
  print(nulldist)
}
save(nulldist,file = 'nulldist_perm_1000.Rdata')


header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
print(header)
#  [1] "Mixture"                      "B cells naive"               
#  [3] "B cells memory"               "Plasma cells"                
#  [5] "T cells CD8"                  "T cells CD4 naive"           
#  [7] "T cells CD4 memory resting"   "T cells CD4 memory activated"
#  [9] "T cells follicular helper"    "T cells regulatory (Tregs)"  
# [11] "T cells gamma delta"          "NK cells resting"            
# [13] "NK cells activated"           "Monocytes"                   
# [15] "Macrophages M0"               "Macrophages M1"              
# [17] "Macrophages M2"               "Dendritic cells resting"     
# [19] "Dendritic cells activated"    "Mast cells resting"          
# [21] "Mast cells activated"         "Eosinophils"                 
# [23] "Neutrophils"                  "P-value"                     
# [25] "Correlation"                  "RMSE"  

load(file = 'nulldist_perm_1000.Rdata')
print(nulldist)
fivenum(print(nulldist)) 
#[1] -0.079400079 -0.009724867  0.019402797  0.056885990  0.604810742
#minimum, lower-hinge, median, upper-hinge, maximum

output <- matrix()
itor <- 1
mix <- dim(Y)[2] ###取mrna中的样本
pval <- 9999
# 表达矩阵的每个样本，都需要计算一下LM22的比例
#iterate through mix
while(itor <= mix){
  
  ##################################
  ## Analyze the first mixed sample
  ##################################
  
  y <- Y[,itor]
  
  #标准化样本数据集
  y <- (y - mean(y)) / sd(y)
  
  #执行SVR核心算法
  result <- CoreAlg(X, y)
  
  #获得结果
  w <- result$w
  mix_r <- result$mix_r
  mix_rmse <- result$mix_rmse
  
  #计算p-value
  if(pval > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
  
  #输出output
  out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
  if(itor == 1) {output <- out}
  else {output <- rbind(output, out)}
  itor <- itor + 1
  
}
head(output)

#save results
write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

#return matrix object containing all results
obj <- rbind(header,output)
obj <- obj[,-1]
obj <- obj[-1,]
obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
rownames(obj) <- colnames(Y)
colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
obj
save(obj,file = 'output_obj.Rdata')

#return matrix object containing all results
obj <- rbind(header,output)
obj <- obj[,-1]
obj <- obj[-1,]
obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
rownames(obj) <- colnames(Y)
colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
obj[1:4,1:4]
#                              B cells naive B cells memory Plasma cells
# TCGA-DK-AA74-01A-11R-A39I-07     0.0125034              0   0.00000000
# TCGA-DK-A3IM-01A-11R-A20F-07     0.0000000              0   0.00000000
# TCGA-GU-A42P-01A-11R-A23W-07     0.1414380              0   0.03075805
# TCGA-4Z-AA7W-01A-11R-A39I-07     0.0000000              0   0.03125399
#                              T cells CD8
# TCGA-DK-AA74-01A-11R-A39I-07  0.11688532
# TCGA-DK-A3IM-01A-11R-A20F-07  0.04647556
# TCGA-GU-A42P-01A-11R-A23W-07  0.02934341
# TCGA-4Z-AA7W-01A-11R-A39I-07  0.37869645

save(obj,file = 'output_obj.Rdata')


load(file = 'input.Rdata')
Y[1:4,1:4]
X[1:4,1:4]
dim(X)
dim(Y)
library(preprocessCore)
library(parallel)
library(e1071)
load(file = 'output_obj.Rdata')

# Step3:将CIBERSORT_Result挑选整理，去除没有差异表达的细胞。
library(dplyr)
library(tidyr)
library(tidyverse)
cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
  rename("Patients" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE"))
# 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Patiens”并去除了后三列。
#并赋值给cibersort_raw。

cibersort_tidy <- cibersort_raw %>%
  remove_rownames() %>%
  column_to_rownames("Patients")
# 将cibersort_raw第一列变为列名后赋值给cibersort_tidy。

flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) < 
                dim(cibersort_tidy)[1]/2)
# 筛选出0值超过样本的一半的一些细胞
cibersort_tidy <- cibersort_tidy[,which(flag)] %>%
  as.matrix() %>%
  t()
# 留下在大部分样本中有所表达的细胞。

bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
# breaks用来定义数值和颜色的对应关系。

# Step4:将CIBERSORT_Result进行可视化
#1）热图
library(pheatmap)
library(RColorBrewer)
pheatmap(
  cibersort_tidy,
  breaks = bk,
  cluster_cols = T,
  scale = "row",
  cluster_row = T,
  border_color = NA,
  show_colnames = F,
  show_rownames = T,
  color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
            colorRampPalette(colors = c("white","red"))(length(bk)/2)
  ))
#调整参数让热图更加美观。

#柱状图可视化细胞占比预测
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))###8种颜色
cibersort_barplot <- cibersort_raw %>%
  gather(key = Cell_type,value = Proportion,2:23)
#使用RColorBrewer包配置需要的色彩方案，使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值给cibersort_barplot

#cibersort_barplot$Patient1 <- factor(cibersort_barplot$Patient,
#                                   levels = str_sort(unique(cibersort_barplot$Patient),
#                                                      numeric = T))

ggplot(cibersort_barplot,aes(Patients,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))
#调整参数让柱状图更加美观。

#直观箱线图
ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,coulour = "black") + theme_bw() + 
  labs(x = "", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypalette(23))
#调整参数让柱状图更加美观。
#这里有个比较新奇的函数是library(RColorBrewer)，是样本配色更加好看

mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
#这里表示按照Set1的配色的前8种
mypalette_1<-brewer.pal(8,"Set1")
image(1:8,1,as.matrix(1:8),col=mypalette_1,xlab="Greens (sequential)",
      ylab="",xaxt="n",yaxt="n",bty="n")

#set1里面配色最多9种，如果设置为10的话

mypalette_1<-brewer.pal(10,"Set1")
image(1:10,1,as.matrix(1:10),col=mypalette_1,xlab="Greens (sequential)",
      ylab="",xaxt="n",yaxt="n",bty="n")
