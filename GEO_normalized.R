## 01.下载GEO数据
rm(list = ls())
library(GEOquery)
f <- "GSE54839.Rdata"
gset = getGEO('GSE54839', destdir=".",getGPL = F)
save(gset,file = f)
## 02.获取ExpressionSet对象，包括的表达矩阵和临床信息
gset=gset[[1]]
exprSet=exprs(gset) ##获取表达矩阵
pdata=pData(gset) ## 获取临床信息
group_list=c(rep('con',30),rep('cocaine',30))
group_list=factor(group_list)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
## 03.校正表达数据
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

## 04.获取平台对应注释R包（https://www.jianshu.com/p/f6906ba703a0）
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
## 05.探针转换以及去重
str(exprSet)
exprSet2 <- exprSet %>% 
  select(-probe_id) %>%  #去掉多余信息
  select(symbol,1:(ncol(exprSet)-1)) %>%  #重新排列，
  mutate(rowMean = rowMeans(exprSet[,1:(ncol(exprSet)-2)])) %>% 
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) #反向选择去除rowMean

save(exprSet2,file = "exprSet2_symbol.Rda")
## 06.差异分析
library(limma)
load("exprSet2_symbol.Rda")
#differential差异分析
##构建分组矩阵
class <- c(rep('con',30),rep('cocaine',30)) 
design <- model.matrix(~factor(class))
colnames(design) <- c("Con","Cocaine")
design

rownames(exprSet2)=exprSet2$symbol
exprSet2=exprSet2[,-1]
#线性模型拟合
fit <- lmFit(exprSet2,design)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=100000) #选择20万是为了输出所有基因
save(allDiff,file = "allDiff.Rda")
#写入表格
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#找出差异两倍以上，pvalue小于0.05，写入表格
diffLab <- subset(allDiff,(logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05)

write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)
#这个是留给热图用的，要不然太长，不好看
diffLab_h <- subset(allDiff,(logFC> 2 | logFC< (-2)) & adj.P.Val < 0.01)

## 07.heatmap热图
#用行名提取数据
library(pheatmap)
#制作一个分组信息
heatdata <- exprSet2[rownames(diffLab),]
#heatdata <- exprSet2[rownames(diffLab_h),]
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=T, # 显示注释
         show_rownames = F,# 显示行名
         scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         border_color = F,
         fontsize = 6)


## 08.volcano火山图
##try drowing a volcano by gglot2
library(ggplot2)
library(ggrepel)
library(dplyr)
data <- allDiff
data$significant <- as.factor(data$P.Value<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("black","red"))+
  labs(title="Volcanoplot", x="log2 (fold change)",y="-log10 (p-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  geom_text_repel(data=subset(data, abs(logFC) > 1), aes(label=gene),col="black",alpha = 0.8)

#ggsave("vocanol.pdf",,width = 7.09, height =5.6,dpi = 300)

## 09.GO分析
## 加载R包
suppressMessages(library(clusterProfiler))
#获得基因列表
gene <- rownames(allDiff)
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de <- gene$ENTREZID
## GO分析
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont="all")
library(ggplot2)
#点图
#pdf(file="dotplot.pdf",width = 10,height = 8)
p1 <- dotplot(go, showCategory = 10,split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p1
#柱状图
#pdf(file="barplot.pdf",width = 10,height = 8)
p2 <- barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
p2
## 10.KEGG富集分析
EGG <- enrichKEGG(gene= gene$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
dotplot(EGG)
