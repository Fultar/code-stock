library(data.table)


rm(list = ls())
gc()

### 数据预处理
file_28_hypo <- fread('../缺氧/28_results/A3SS.MATS.JC.FDR0.05.txt')
colnames(file_28_hypo)
data_hypo <- file_28_hypo[,c(2,6,7)]
colnames(data_hypo) <- c('GeneID','Start','End')
as_name <- paste0(data_hypo$GeneID,'_',data_hypo$Start,'_',data_hypo$End)
data_counts <- data.frame(GeneID = as_name, Control_counts = file_28_hypo$IJC_SAMPLE_1,
                          Treat_counts = file_28_hypo$IJC_SAMPLE_2)


### 导入gtf文件
##导入数据与计算
library(GenomicFeatures) 
library(tidyverse)
txdb <- makeTxDbFromGFF("../database/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf",format="gtf") # 将数据导入为TxDb对象
exons_gene <- exonsBy(txdb, by = "gene") 
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))}) #计算总外显子长度
exons_gene_lens[1:10]

##转换为dataframe
gene_length <- sapply(exons_gene_lens,function(x){x})
id_length <- as.data.frame(gene_length)
id_length

## 到ensembl BioMart下载geneid和genename对应的文件
## biomart > 找到罗非鱼物种名 > 只把geneid 和genename打勾 > TSV格式导出 > 把表头的空格删除

GN <- fread('tilapia_mart_export.txt',header = T)

# 合并数据
rownames(GN) <- GN$GenestableID #用geneid作为行名
id_GN_length <- cbind(GN,id_length) #合并一下
colnames(id_GN_length) <- c("gene_id","gene_name","length") #改一下行名
id_GN_length #得到最终的dataframe

write.table(id_GN_length,'Tilapia_gene_length.txt',row.names = F, sep = '\t', quote = F)



### 表达量计算
library(dplyr)
length <- distinct(length, gene_name, .keep_all = T)
mergecount <- merge(rawcount, length, by = 'gene_name')
FPKMlength <- mergecount[,c(1, ncol(mergecount))]
FPKMcount <- mergecount[, -c(2, ncol(mergecount))]
rownames(FPKMcount) <- FPKMcount$gene_name
FPKMcount <- FPKMcount[,-1]

## 计算 TPM ##

kb <- FPKMcount$Length / 1000
kb
countdata <- FPKMcount[,5:7]   #r的索引是从1开始的，5：7选择的是count里面每个样本对应的reads数的列
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
#将上面计算好的tpm保存到本地
ootpm <- as.data.frame(tpm)
write.csv(ootpm, file="C:/Users/Desktop/ootpm.csv",quote=FALSE)


## 计算 FPKM ##
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
ofpkm <- as.data.frame(fpkm)
write.csv(ofpkm , file="C:/Users/Desktop/2tpm.csv",quote=FALSE)