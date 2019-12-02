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

##获取平台对应注释R包（https://www.jianshu.com/p/f6906ba703a0）
library("illuminaHumanv3.db")
#获取探针
probe2symbol_df <- toTable(get("illuminaHumanv3SYMBOL"))
#使用相同的名称
names(probe2symbol_df)[1]
exprSet=data.frame(exprSet,check.names = F)
exprSet$probe_id=row.names(exprSet)
library(dplyr)
exprSet=dplyr::inner_join(exprSet,probe2symbol_df,by="probe_id")
#看一下symbol有没有重复，发现只有41932个，所以需要去重
###探针转换以及去重
str(exprSet)
exprSet2 <- exprSet %>% 
  select(-probe_id) %>%  #去掉多余信息
  select(symbol,1:(ncol(exprSet)-1)) %>%  #重新排列，
  mutate(rowMean = rowMeans(exprSet[,1:(ncol(exprSet)-2)])) %>% 
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) #反向选择去除rowMean

save(exprSet2,file = "exprSet2_symbol.Rda")
