library(ggplot2)
library(data.table)

rm(list = ls())
gc()

### 提取所需基因名
vene_inter <- fread('all_venn_inter.txt')

co_genes <- vene_inter[1,2]
co_genes <- c(t(co_genes))
co_genes <- strsplit(co_genes,', ')[[1]]
co_genes


### 合并文件
{file_28hypo <- fread('../缺氧/28_results/ALL_DAS.MATS.JC.txt')
file_28hypo$Group <- c('Hypoxia_28')

file_32hypo <- fread('../缺氧/32_results/ALL_DAS.MATS.JC.txt')
file_32hypo$Group <- c('Hypoxia_32')

file_36hypo <- fread('../缺氧/36_results/ALL_DAS.MATS.JC.txt')
file_36hypo$Group <- c('Hypoxia_36')

file_28vs32 <- fread('../温度/28VS32/ALL_DAS.MATS.JC.txt')
file_28vs32$Group <- c('Temperature_28vs32')

file_28vs36 <- fread('../温度/28VS36/ALL_DAS.MATS.JC.txt')
file_28vs36$Group <- c('Temperature_28vs36')

file_32vs36 <- fread('../温度/32VS36/ALL_DAS.MATS.JC.txt')
file_32vs36$Group <- c('Temperature_32vs36')


file_merge <- rbind(file_28hypo,file_32hypo,file_36hypo,
                    file_28vs32,file_28vs36,file_32vs36)
}


### 提取含目标基因的行
as_list <- list()
for (gene in co_genes){
  as_list[[gene]] <- file_merge[grep(gene,file_merge$GeneID),]
}
library(writexl)
writexl::write_xlsx(as_list, '缺氧和高温组共有的21个可变剪切基因的数据统计.xlsx')
