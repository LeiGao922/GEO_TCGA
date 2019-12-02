#下载GEO数据
rm(list = ls())
library(GEOquery)
f <- "GSE54839.Rdata"
gset = getGEO('GSE54839', destdir=".",getGPL = F)
save(gset,file = f)
## 获取ExpressionSet对象，包括的表达矩阵和临床信息
gset=gset[[1]]
exprSet=exprs(gset) ##获取表达矩阵
pdata=pData(gset) ## 获取临床信息
group_list=c(rep('con',30),rep('cocaine',30))
group_list=factor(group_list)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
## 校正表达数据
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
## 判断是否需数据校正
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

