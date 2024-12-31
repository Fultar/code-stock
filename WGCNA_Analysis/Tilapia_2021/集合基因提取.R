library(data.table)

# 温度组
file_temp_venn <- fread('temperature_venn_inter.txt')

gene_28vs36 <- file_temp_venn[4,2]
gene_28vs36 <- data.frame(strsplit(gene_28vs36$GeneID,', '))
colnames(gene_28vs36) <- c('GeneID')


# 温度和缺氧组联合
file_all_venn <- fread('all_venn_inter.txt')
gene_28vs36_hypo <- file_all_venn[32,2]
gene_28vs36_hypo <- data.frame(strsplit(gene_28vs36_hypo$GeneID,', '))
colnames(gene_28vs36_hypo) <- c('GeneID')


## 读取可变剪切事件
file_28vs36_as <- fread('DAS_28VS36.txt')

data_28vs36_hypo_as <- merge(file_28vs36_as, gene_28vs36_hypo,by = 'GeneID')
data_28vs36_as <- merge(file_28vs36_as, gene_28vs36,by = 'GeneID')

# 提取特定AS类别

data_28vs36_as_se <- data_28vs36_as[data_28vs36_as$type == 'SE',]
data_28vs36_as_se <- data_28vs36_as_se[!duplicated(data_28vs36_as_se$GeneID)]
write.table(data_28vs36_as_se, '温度组28度VS36度特有的SE可变剪切基因.txt', row.names = F, quote = F)


data_28vs36_hypo_as_se <- data_28vs36_hypo_as[data_28vs36_hypo_as$type == 'SE',]
data_28vs36_hypo_as_se <- data_28vs36_hypo_as_se[!duplicated(data_28vs36_hypo_as_se$GeneID)]
write.table(data_28vs36_hypo_as_se, '温度组与缺氧组比较-28度VS36度特有的SE可变剪切基因.txt', row.names = F, quote = F)
