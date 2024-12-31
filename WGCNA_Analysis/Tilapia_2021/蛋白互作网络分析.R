library(STRINGdb)

rm(list = ls())
gc()


# 初始化 STRINGdb，指定物种（9606 是人类的 taxonomy ID）
string_db <- STRINGdb$new(version = "11.5", species = 8128, 
                          score_threshold = 100, input_directory = "")

# 获取 STRINGdb 数据
# 输入一个蛋白质的基因符号列表，或者可以获取全局的 PPI 网络
file_gene <- fread('DAS_36_hypoxia.txt')
genes <- file_gene[c(1:10),1]
ppi_data <- string_db$get_interactions(genes)

# 查看 PPI 数据
head(ppi_data)
