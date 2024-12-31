# 魔幻操作，一键清空
rm(list = ls())


library(stringr) # 标签换行
library(AnnotationDbi)
library(DOSE)
library(ggplot2) # 绘图
library(ggrepel) # 标签相关
library(ggsci) #配色包

# 安装字体
library(showtext)
showtext_auto(enable = TRUE)
font_add('SimSun', 'simsun.ttc')

##对于有参考基因组物种的分析，可以在相关软件包中直接加载该物种的背景基因集
#（1）对于常见的模式物种，例如人类，有些专门的 R 包数据库，例如人类参考基因组 hg19 的
library(org.Hs.eg.db)

#（2）对于不常见的物种，但却是存在参考基因组的情况
#通过 AnnotationHub 包索引基因组，例如我先前分析过一个绵羊（Ovis aries）的，这样导入基因数据库
#library(AnnotationHub)
#hub <- AnnotationHub()
#query(hub, 'Ovis aries')  #输入绵羊（Ovis aries）的名称进行匹配
#sheep <- hub[['AH72269']]  #返回了数据库编号 AH72269，就可以加载该库

##clusterProfiler的 GO 富集
library(clusterProfiler)

#读取基因列表文件中的基因名称
#genes <- read.csv('up-gene.csv', header = TRUE)
#id_list <- genes[,1]
# 将gene symbol转换为Entrez ID,防止分析出错
#id_list <- mapIds(org.Hs.eg.db,id_list,"ENTREZID","SYMBOL")
# 去除转换不成功的结果，即id=NA的情况
#id_list <- na.omit(id_list)
#GO富集分析
#enrich.go <- enrichGO(gene = id_list,  #待富集的基因列表
                      #OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库
                      #keyType = 'ENTREZID',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      #ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      #pAdjustMethod = 'fdr',  #指定 p 值校正方法
                     # pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                     # qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                     # readable = T)

#输出
#write.table(enrich.go, 'enrich.go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
enrich.go <- read.csv('Regulated_GO_enrich_all.csv',header = T)
go.res <- data.frame(enrich.go)
ONTOLOGY <- go.res[,1]
goBP <- subset(go.res,subset = (ONTOLOGY == "Biological Process"))[1:10,] 
goCC <- subset(go.res,subset = (ONTOLOGY == "Cellular Component"))[1:10,] 
goMF <- subset(go.res,subset = (ONTOLOGY == "Molecular Function"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)%>% na.omit(goMF)
# 使画出的GO term的顺序与输入一致
go.df$GO.Terms.Description <- factor(go.df$GO.Terms.Description,
                                     levels = rev(go.df$GO.Terms.Description))
go.df$NO <- factor(go.df$NO,levels = c(30:1))
### 绘图
##条形图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = NO, 
                     y = X.LOG10.P.Value.,
                     fill = GO.Terms.Level.1))+ # 横轴坐标及颜色分类填充
  #geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  geom_col(width = 0.6,
           show.legend = FALSE)+
  #scale_fill_manual(values = c('#C0392B', '#13EF4F', '#1667F3'))+
  scale_color_aaas()+
  scale_fill_aaas()+ #使用配色方案
  coord_flip()+
  theme_bw()+ # 横纵坐标反转及去除背景色
  #scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "-lgP value",title = "GO Enrichment BarPlot")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), # 图边距
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  facet_grid(GO.Terms.Level.1~., scale = 'free_y', space = 'free_y')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
go_bar
ggsave(go_bar,filename = "GO_Barplot_Numbers.pdf",width = 6,height = 7)



### KEGG富集分析
#1、KEGG富集
#kk <- enrichKEGG(gene = id_list,keyType = "kegg",organism= "human", qvalueCutoff = 0.5, pvalueCutoff=0.5)

#2、可视化
###柱状图
#保存KEGG结果
library(DOSE)
#kegg_results<-setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#hh <- as.data.frame(kegg_results)#自己记得保存结果哈！
#write.table(hh,"KEGG_results.txt",sep = '\t', row.names = FALSE, quote = FALSE)
rm(list = ls())

hh <- read.csv('Regulated_kegg_enrich_all.csv',header = T)
rownames(hh) <- 1:nrow(hh)
hh$Fold.Enrichment <- sort(hh$Fold.Enrichment,decreasing = T)

##KEGG条形图
file_tilapia_class <- read.table('Tilapia_Regulated-KEGG_map_classify.txt',header = T,sep = '\t')
file_carp_class <- read.table('GrassCarp_Regulated-KEGG_map_classify.txt',header = T,sep = '\t')

data_tilapia_class <- file_tilapia_class
data_tilapia_class$Num <- c(1:31)
data_tilapia_class$Num <- factor(data_tilapia_class$Num, levels = c(1:31))

kegg_bar <- ggplot(data_tilapia_class,
                   aes(y = Number.of.proteins,
                      x = Num,
                      fill = KEGG.map.level_1))+
  geom_bar(stat = "identity")+  ####柱子宽度
  coord_flip()+  ##颠倒横纵轴
  #scale_fill_gradient()+ #颜色自己可以换
  labs(title = "KEGG Map Classify",
       x = "Map", 
       y = "Protein Number")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),# 去除网格线
        axis.line.y = element_line(size = 1, color = "black"),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.title.x = element_text(size=18,face = "bold",color = "black"), # 调整坐标轴字体
        axis.title.y = element_text(size=18,face = "bold",color = "black"),
        axis.text = element_text(size = 14,color = "black"),
        axis.text.x = element_text(hjust = 0.5), # 调整x轴刻度
        axis.line = element_line(size = 1.5),
        axis.ticks.y = element_line(color="black",size=1,lineend = 10),
        axis.ticks.x = element_line(color="black",size=1,lineend = 10),
        panel.border = element_blank(),
        legend.key.size = unit(15, "pt"))+
  scale_y_discrete(expand = c(0,0), breaks = c(rep(0:120),25))
kegg_bar
ggsave(kegg_bar,filename="kegg_bar.pdf",width = 8,height = 5)



###KEGG气泡图
#rownames(hh) <- 1:nrow(hh)
#hh$Fold.Enrichment <- sort(hh$Fold.Enrichment,decreasing = T)
kegg_hh <- hh[c(1:20),]
colnames(kegg_hh)[7] <- c("pvalue")
colnames(kegg_hh)[2] <- c('GeneNumber')
kegg_hh <- arrange(kegg_hh,desc(kegg_hh[,6]))
kegg_hh <- kegg_hh[c(1:15),]
kegg_hh$KEGG.pathway <- factor(kegg_hh$KEGG.pathway,levels = rev(kegg_hh$KEGG.pathway))
kegg_hh$NO <- c(1:15)
kegg_hh$NO <- factor(kegg_hh$NO,levels = c(15:1))

kegg_bubble <- ggplot(kegg_hh,
                      aes(y=NO,
                          x=Fold.Enrichment,
                          colour=-1*log10(pvalue),
                          size=GeneNumber))+
                geom_point()+
                scale_size(range=c(2, 6))+
                scale_color_gradient(high = "#EF204C",
                                     low = "#5D7FF3")+
                labs(color=expression(-log[10](PValue),
                                      size="GeneNumber"), 
                      x="Fold Enrichment",
                      y="Pathways",
                      title="KEGG Pathway Enrichment")+
                theme_bw()+
                theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18),
                      #panel.grid.major = element_blank(),panel.grid.minor = element_blank(),# 去除网格线
                      axis.line.y = element_line(size = 1, color = "black"),
                      axis.line.x = element_line(size = 1, color = "black"),
                      axis.title.x = element_text(size=18,face = "bold",color = "black"), # 调整坐标轴字体
                      axis.title.y = element_text(size=18,face = "bold",color = "black"),
                      axis.text = element_text(size = 18,face = "bold",color = "black"),
                      axis.text.x = element_text(hjust = 0.5), # 调整x轴刻度
                      axis.line = element_line(size = 1.5),
                      axis.ticks.y = element_line(color="black",size=1,lineend = 10),
                      axis.ticks.x = element_line(color="black",size=1,lineend = 10),
                      #panel.border = element_blank(),
                      legend.key.size = unit(15, "pt"))+
  scale_x_continuous(expand = c(0, 0),limits = c(2,7),breaks = c(1:7))
                      #text = element_text(family = 'SimSun'))
kegg_bubble
ggsave(kegg_bubble,filename="kegg_bubble_Numbers.pdf",width = 6,height = 6)
write.table(kegg_hh,file = 'kegg_bubble.txt')
