library(reshape2)
library(ggplot2)

rm(list = ls())
gc()


file_wgcna <- read.table('Module_gene_name_greenyellow_temperature.txt', header = T)
colnames(file_wgcna)[1] <- c('GeneID')

file_as <- read.table('../数据分析/DAS_28VS36.txt',header = T)

data_cogene <- merge(file_as,file_wgcna,by = 'GeneID')

file_tmp <- read.csv('tilapia_gene_tpm_log2+1.csv',header = T)
colnames(file_tmp)[1] <- c('GeneID')

data_cogene <- merge(data_cogene, file_tmp, by = 'GeneID')

data_cogene2 <- melt(data_cogene)

ggplot(data_cogene2, aes(x = variable, y = GeneID, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 5) +  # 使用蓝-白-红的渐变
  labs(x = "Columns", y = "Rows", fill = "Value") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 25),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),# 去除网格线
        axis.line.y = element_line(size = 1, color = "black"),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.title.x = element_text(size = 18,face = "bold",color = "black"), # 调整坐标轴字体
        axis.title.y = element_text(size = 18,face = "bold",color = "black"),
        axis.text = element_text(size = 14,face = "bold",color = "black"),
        axis.text.x = element_text(hjust = 0.5, angle = 90), # 调整x轴刻度
        
        axis.line = element_line(size = 1.5),
        axis.ticks.y = element_line(color = "black",size=1,lineend = 10),
        axis.ticks.x = element_line(color = "black",size=1,lineend = 10),
        #legend.position = "none",
        legend.key.size = unit(20, "pt"), # 图例大小
        panel.border = element_blank())
ggsave('温度性状-WGCNA黄绿色模块-28VS36可变剪切-共有基因的表达量热图.pdf',
       width = 18, height = 12)
write.table(data_cogene,'温度性状-WGCNA黄绿色模块-28VS36可变剪切-共有基因及表达量.txt', quote = F)

###---------------缺氧性状可变剪切与WGCNA模块---------------####
library(data.table)

### darkgrey module
file_28hy <- fread('../可变剪切分析/DAS_28_results.txt')
file_28hy$Group <- c('Hypoxia_28')
file_28hy <- file_28hy[,c(2,4,5,14,17,18,19)]

file_32hy <- fread('../可变剪切分析/DAS_32_results.txt')
file_32hy$Group <- c('Hypoxia_32')
file_32hy <- file_32hy[,c(2,4,5,14,17,18,19)]

file_36hy <- fread('../可变剪切分析/DAS_36_results.txt')
file_36hy$Group <- c('Hypoxia_36')
file_36hy <- file_36hy[,c(2,4,5,14,17,18,19)]


data_all_hy <- rbind(file_28hy, file_32hy, file_36hy)

file_grey <- fread('Module_gene_name_darkgrey_temperature.txt')
colnames(file_grey) <- c('GeneID')
data_cogene <- merge(file_grey, data_all_hy, by = 'GeneID')
#data_cogene <- data_cogene[!duplicated(data_cogene$GeneID),]


write.table(data_cogene, file = '缺氧性状可变剪切与黑灰色模块共有基因.txt', row.names = F, quote = F, sep = '\t')



### magenta module
file_magenta <- fread('Module_gene_name_magenta_hypoxia.txt')
colnames(file_magenta) <- c('GeneID')
data_cogene <- merge(file_magenta, data_all_hy, by = 'GeneID')
#data_cogene <- data_cogene[!duplicated(data_cogene$GeneID),]

write.table(data_cogene, file = '缺氧性状可变剪切与品红色模块共有基因.txt', row.names = F, quote = F,sep = '\t')

###---------------耐热性状可变剪切与WGCNA模块---------------####
rm(list = ls())
gc()

file_28vs36 <- fread('../可变剪切分析/DAS_28VS36.txt')
file_28vs36$Group <- c('Temperature_28vs36')
file_28vs36 <- file_28vs36[,c(2,4,5,14,17,18,19)]


file_28vs32 <- fread('../可变剪切分析/DAS_28VS32.txt')
file_28vs32$Group <- c('Temperature_28vs32')
file_28vs32 <- file_28vs32[,c(2,4,5,14,17,18,19)]


file_32vs36 <- fread('../可变剪切分析/DAS_32VS36.txt')
file_32vs36$Group <- c('Temperature_32vs36')
file_32vs36 <- file_32vs36[,c(2,4,5,14,17,18,19)]


data_all_temp <- rbind(file_28vs32, file_28vs36, file_32vs36)
#data_all_temp <- data_all_temp[!duplicated(data_all_temp),]

### greenyellow module
file_greenyellow <- fread('Module_gene_name_greenyellow_temperature.txt')
colnames(file_greenyellow) <- c('GeneID')

data_cogene <- merge(data_all_temp, file_greenyellow, by = 'GeneID')

write.table(data_cogene, file = '耐热性状可变剪切与黄绿色模块共有基因.txt', row.names = F, quote = F,sep = '\t')

### darkred
file_darkred <- fread('Module_gene_name_darkred_temperature.txt')
colnames(file_darkred) <- c('GeneID')

data_cogene <- merge(data_all_temp, file_darkred, by = 'GeneID')

write.table(data_cogene, file = '耐热性状可变剪切与深红色模块共有基因.txt', row.names = F, quote = F,sep = '\t')





