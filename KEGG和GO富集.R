rm(list = ls())
gc

library(data.table)

### 1.数据整理
{file_36_A3SS <- fread('../36_results/A3SS.MATS.JC.FDR0.05.txt', header=T)
file_36_A5SS <- fread('../36_results/A5SS.MATS.JC.FDR0.05.txt', header=T)
file_36_RI <- fread('../36_results/RI.MATS.JC.FDR0.05.txt', header=T)
file_36_MXE <- fread('../36_results/MXE.MATS.JC.FDR0.05.txt', header=T)
file_36_SE <- fread('../36_results/SE.MATS.JC.FDR0.05.txt', header=T)

file_36_A3SS$Type <- c('A3SS')
file_36_A5SS$Type <- c('A5SS')
file_36_RI$Type <- c('RI')
file_36_MXE$Type <- c('MXE')
file_36_SE$Type <- c('SE')

data_36_A3SS <- file_36_A3SS[,c(2,24)]
data_36_A5SS <- file_36_A5SS[,c(2,24)]
data_36_MXE <- file_36_MXE[,c(2,26)]
data_36_RI <- file_36_RI[,c(2,24)]
data_36_SE <- file_36_SE[,c(2,24)]

data_36_AS <- rbind(data_36_A3SS,data_36_A5SS,data_36_MXE,data_36_RI,data_36_SE)

write.table(data_36_AS,'ALL_DAS_36.txt',sep = '\t', quote = F, row.names = F)}

### 2.GO富集
rm(list = ls())
gc()

library(AnnotationDbi)
#library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(data.table)

# 查询罗非鱼的物种注释包
library(AnnotationHub)
hub <- AnnotationHub()
#query(hub, 'Oreochromis niloticus')  
tilapia <- hub[['AH114368']]
length(keys(tilapia))
columns(tilapia)

###### GO富集
module <- 'magenta'
pheno <- 'hypoxia'
{#module <- 'darkred'
#pheno <- 'temperature'
filename <- c(paste0('Module_gene_name_',module,'_', pheno))  ## 选择要分析的组别
data_AS <- fread(paste0(filename, ".txt"))
#id_list <- data_AS$GeneID
id_list <- data_AS$x

# 转换ID
trans_id <- bitr(id_list,'SYMBOL', c("ENTREZID","REFSEQ", "GO", "ONTOLOGY"), tilapia)

# 去除转换不成功的结果，即id=NA的情况
trans_id <- na.omit(trans_id)

# GO富集
ont_name <- c('ALL')   #选择生物通路大类：CC/BP/MF/ALL

enrich.go <- enrichGO(gene = trans_id[,2],  
                      OrgDb = tilapia, 
                      keyType = 'ENTREZID',  
                      ont = ont_name,          #选择生物通路大类：CC/BP/MF/ALL
                      pAdjustMethod = 'BH',  
                      pvalueCutoff = 1,  
                      qvalueCutoff = 1,  
                      readable = T)
# 需要选取单通路类型才能运行这行代码
#clusterProfiler::goplot(enrich.go)

#ggsave(paste0(filename,'_enrichGO_',ont_name,'.png'), height = 12, width = 16)

write.table(enrich.go, paste0(filename, '_enrichGO.txt'), 
            sep = '\t', row.names = FALSE, quote = FALSE)
}


#### GO富集绘图
{library(ggsci)
#enrich.go <- read.table('DAS_28VS36_enrichGO.txt', header = T, sep = '\t')
go.res <- enrich.go
ONTOLOGY <- go.res[,1]
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:10,] 
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:10,] 
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
#go.df$NO <- c(1:nrow(go.df))
#go.df$NO <- factor(go.df$NO,levels = c(nrow(go.df):1))
### 绘图
##条形图
#go.df[c(29),c(3)] <- c(paste0('oxidoreductase activity,','tacting on paired donors,with incorporation or reduction of molecular oxygen,','\n', 'reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen'))
colors <- c('#3498db', '#e67e22', '#db3457')

go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = reorder(Description, -pvalue), 
                     y = -log(pvalue),
                     fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  #geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  geom_col(width = 0.6,
           show.legend = FALSE)+
  scale_fill_manual(values = colors)+
  #scale_color_aaas()+
  #scale_fill_aaas()+ #使用配色方案
  coord_flip()+
  theme_bw()+ # 横纵坐标反转及去除背景色
  #scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "-lgP value",title = "GO Enrichment BarPlot")+ # 设置坐标轴标题及标题
  # 条目名在柱状图里
  geom_text(aes(y = 0.05, label = Description,), size = 4, hjust = 0) +
  # 显示条目基因
  #geom_text(aes(y = 0.1, label = geneID, ), size = 3.5,fontface = 'italic', hjust = 0, vjust = 1.8) +
  theme(
        axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(color = 'black', size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), # 图边距
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        # 隐藏y轴标签
        axis.text.y = element_blank(),
        # x轴刻度线朝里
        axis.ticks.length.y = unit(-0.2, "cm"),
        axis.ticks.length.x = unit(-0.2, "cm")
        )+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
go_bar

ggsave(paste0(filename, '_enrichGO_TOP30_text.pdf'), width = 16,height = 12)
ggsave(paste0(filename, '_enrichGO_TOP30_text.png'), width = 16,height = 12, dpi = 600)
}

##### KEGG通路富集分析
{
#ilename <- c('DAS_32VS36')  ## 选择要分析的组别
#data_AS <- fread(paste0(filename, ".txt"))
#id_list <- data_AS$GeneID

# 转换ID
#trans_id <- bitr(id_list,'SYMBOL', c("ENTREZID","REFSEQ", "GO", "ONTOLOGY"), tilapia)

# 去除转换不成功的结果，即id=NA的情况
#trans_id <- na.omit(trans_id)

# 分析

kk <- enrichKEGG(gene = trans_id[,2], keyType = "kegg",organism = 'onl', qvalueCutoff = 1, pvalueCutoff= 1,
                  minGSSize = 10, maxGSSize = 5000)
# ID转换
kk2 <- setReadable(kk, OrgDb = tilapia, keyType = "ENTREZID")

write.table(kk2, paste0(filename, '_enrichKEGG.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}
### KEGG通路图
{library('pathview')
pathname <- c('E:/课业文件/研究生/周艳罗非鱼可变剪切/周艳罗非鱼高温缺氧可变剪切/温度/28VS36/')
filename <- c('ALL_DAS.MATS.JC.txt')
file_das <- fread(paste0(pathname,filename))
colnames(file_das)
data_das_pvalue <- file_das[,c(1,2,13)]
data_das_pvalue <- unique(data_das_pvalue, by = 'ID')

geneid <- data_das_pvalue$ID
rownames(data_das_pvalue) <- geneid
head(data_das_pvalue)
data_das_pvalue$log2Pvalue <- -log2(data_das_pvalue$PValue)
data_das_pvalue$expr <- c(0.8)

#### 列名为基因的ID，每列的样本数据可以是表达量、显著性p值等需要区分的数值
onlpathway <- pathview(gene.data  = data[,1],
                     pathway.id = "03040",
                     species    = "onl",
                     limit      = list(gene=max(abs(data[,1])), cpd=1))


}
### KEGG柱状图
{#file_kegg <- fread('36Group_ALL_DAS_enrich_KEGG.txt')
data_kegg <- kk2[c(1:10),]
data_kegg$NO <- c(1:10)

write.table(data_kegg[,-ncol(data_kegg)], paste0(filename, '_enrichKEGG_TOP10.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}

###KEGG气泡图
{library(stringr)
#data_kegg <- fread('36Group_TOP3_enrich_KEGG.txt')
data_kegg <- kk2[c(1:10),]
data_kegg$Description <- str_sub(data_kegg$Description, 1,-40)
kegg_bubble <- ggplot(data_kegg,
                      aes(y = reorder(Description, FoldEnrichment),
                          x = FoldEnrichment,
                          colour = -log(pvalue),
                          size = Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_color_gradient(high = "#db3457",
                       low = "#3498db")+
  labs(color=expression(-log[10](pvalue),
                        size="Count"), 
       x="FoldEnrichment",
       y="Pathways",
       title="KEGG Pathway Enrichment")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18),
        #panel.grid.major = element_blank(),panel.grid.minor = element_blank(),# 去除网格线
        axis.line.y = element_line(size = 1, color = "black"),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.title.x = element_text(size=18,face = "bold",color = "black"), # 调整坐标轴字体
        axis.title.y = element_text(size=18,face = "bold",color = "black"),
        axis.text = element_text(size = 14,face = "bold",color = "black"),
        axis.text.x = element_text(hjust = 0.5), # 调整x轴刻度
        axis.line = element_line(size = 1.5),
        axis.ticks.y = element_line(color="black",size=1,lineend = 10),
        axis.ticks.x = element_line(color="black",size=1,lineend = 10),
        #panel.border = element_blank(),
        legend.key.size = unit(15, "pt"))+
  scale_x_continuous(expand = c(0, 0),limits = c(1,7),breaks = c(1,seq(0,20,1)))
kegg_bubble

ggsave(paste0(filename, '_enrichKEGG_Bubble.pdf'),width = 12,height = 5)
ggsave(paste0(filename, '_enrichKEGG_Bubble.png'),width = 12,height = 5, dpi = 500)
}






########### ------------------测试数据------------------ ##########


library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)
library(org.Hs.eg.db)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

goplot(ego)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

data(gse16873.d)
head(gse16873.d)
p <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa",
              out.suffix = "gse16873", kegg.native = T)
