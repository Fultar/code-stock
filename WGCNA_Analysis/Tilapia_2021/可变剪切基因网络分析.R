rm(list = ls())
gc()

BiocManager::install('mGSZ')

library(ASpediaFI)

# 前期准备文件：bam经可变剪切分析后的剪切事件，这里我用rMATS计算好了，ASpediaFI包也可以计算剪切事件
# 可变剪切文件需要PSI值，即包含该剪切事件的mRNA比例。rMATS分析文件中IncLevel1和IncLevel2为两个样本组的PSI值。
# 此外，还需要样本的表达量文件，可用limma或DESeq2等包分析基因表达量





















