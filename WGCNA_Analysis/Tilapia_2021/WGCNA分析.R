library(WGCNA)
library(data.table)
library(dbplyr)

rm(list = ls())
gc()

#############-------------------根据DESeq2生产的counts计算TPM--------------############
file_expr <- fread('all_des_output.csv')
colnames(file_expr)

file_feature <- fread('all_feature.txt')

# 获取需要的列：第1列-基因名；第2-第n列-样品的counts数；最后一列-每个基因的长度
data_counts <- file_expr[,c(-2:-7)]
colnames(data_counts)[1] <- c('Geneid')
data_counts <- na.omit(data_counts)
# 删去表达量为0的行
row_sums <- rowSums(data_counts == 0)

data_counts <- data_counts[row_sums <= 3, ]

data_length <- file_feature[,c(1,6)]
# 合并
data_counts_filtered <- data.frame(merge(data_counts, data_length, by = 'Geneid'))
row.names(data_counts_filtered) <- data_counts_filtered$Geneid
data_counts_filtered <- data_counts_filtered[,-1]
head(data_counts_filtered)

#计算tpm
## 前处理
head(data_counts_filtered)

kb <- data_counts_filtered$Length / 1000
head(kb)

# 计算TPM 
# 【至少要有两列Row reads count数据，才能计算，否则报错。
# 其中的mycounts[,1:2]根据你的样品数量，像这里我只有两个样品就是1：2，如果你是10个样品就是 1：10】
countdata <- data_counts_filtered[,1:18]

head(countdata)

rpk <- countdata / kb
head(rpk)

tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)

tpm_log <- log2(tpm + 1)

write.csv(tpm, 'tilapia_gene_tpm.csv', quote = F)
write.csv(tpm_log, 'tilapia_gene_tpm_log2+1.csv', quote = F)



#############-------------------WGCNA分析--------------############
rm(list = ls())
gc()


load('networkConstruction.RData')
load('WGCNA.RData')

# 表达矩阵：常规表达矩阵即可，数据可以先经log2(TPM+1)处理均一化
exp_tpm <- read.csv("tilapia_gene_tpm_log2+1.csv", row.names = 1)
exp_tpm <- t(exp_tpm)

head(exp_tpm)[1:5, 1:5]

# 判断一下数据质量
# 通过goodSamplesGenes函数生成一个质量检查结果对象gsg
gsg <- goodSamplesGenes(exp_tpm, verbose = 3)

# `verbose` 参数被设置为 `3`，用于控制函数 `goodSamplesGenes` 的详细输出程度。
# `verbose = 0`：不产生任何输出，只返回结果，通常用于静默模式。
# `verbose = 1`：产生基本的信息输出，以提供一些关于函数执行进度的信息。
# `verbose = 2`：产生更详细的输出，可能包括一些中间步骤的信息。
# `verbose = 3`：产生最详细的输出，通常包括每个步骤的详细信息，用于调试和详细分析。

# 检查质量检查结果对象中的allOK属性，如果为FALSE，表示质量不好
gsg$allOK

# 如果质检结果为FALSE，运行以下代码过滤：
if (!gsg$allOK) 
{ 
  # 如果goodGenes属性中有不好的基因，打印并移除它们
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(exp_tpm)[!gsg$goodGenes], collapse = ", "))) 
  
  # 如果goodSamples属性中有不好的样本，打印并移除它们
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(exp_tpm)[!gsg$goodSamples], collapse = ", "))) 
  
  # 根据质量检查结果对象中的goodSamples和goodGenes属性，过滤数据集exp_tpm
  exp_tpm <- exp_tpm[gsg$goodSamples, gsg$goodGenes] 
}

# 第一步，咱们就是构建样本的系统聚类树
# 主要为了查看是否有离群样本
sampleTree <- hclust(dist(exp_tpm), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
# 这里样本较少我没有过滤离群样本

# 绘制样本聚类图（上）与样本性状热图（下）

# 这一步是为了将数值转换为颜色编码，用于可视化，将数值映射到颜色，我们可以更直观地展示数据的特征或模式
drug_screen_parp <- read.csv('phenotype.csv',row.names = 1)
traitColors <- numbers2colors(drug_screen_parp$temperature, signed = FALSE)
head(traitColors)

# 可视化
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(drug_screen_parp),
                    main = "Sample dendrogram and trait heatmap")

# 我们可以保存一下处理好的表达矩阵和临床信息（也就是所谓的性状矩阵）。
datExpr <- exp_tpm %>% filter(rownames(exp_tpm) %in% rownames(drug_screen_parp))
dim(datExpr)
saveRDS(datExpr, file = "datExpr.rds")

# 挑选最佳软阈值

# 设置 power 参数选择范围，可以自行修改设定
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sft <- pickSoftThreshold(exp_tpm, powerVector = powers, verbose = 5)

# 我们看一下软阈值选择结果
sft

# 绘制软阈值和拟合指标的关系图
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 1, col = "red")

# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.85, col = "red")



# 绘制软阈值对平均连接度的影响图
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 1, col = "red")



# 正式构建加权共表达网络

# 一步法，这步时间可能有点久
nGenes = ncol(exp_tpm)
nSamples = nrow(exp_tpm)
net <- blockwiseModules(exp_tpm, power = sft$powerEstimate, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)
# power = sft$powerEstimate：这里指定了软阈值的取值。sft$powerEstimate 就是我们之前计算得到的软阈值，它用来控制共表达网络的连接强度。
# maxBlockSize = nGenes：指定了最大的模块大小。在构建基因模块时，会将基因分成多个子模块，这个参数用来限制子模块的最大大小。nGenes 是基因的总数，通常用来作为最大模块大小的上限。（4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个，计算资源允许的情况下最好放在一个block里面。）
# TOMType = "unsigned"：指定了共表达网络的类型，这里设置为 "unsigned"，表示构建无符号的网络。计算TOM矩阵时，是否考虑正负相关性，默认为"signed"，但是根据幂律转换的邻接矩阵(权重)的非负性，所以我们认为这里选择"signed"也没有太多的意义。
# minModuleSize = 30：指定了每个模块的最小大小，也就是每个基因模块中至少包含多少个基因。
# reassignThreshold = 0：控制模块的重新分配。如果模块中的基因与其他模块的相关性低于这个阈值，它们可能会被重新分配到其他模块。
# mergeCutHeight = 0.25：用来合并基因模块。当两个模块的相关性超过这个阈值时，它们可能会被合并成一个更大的模块。越大得到的模块就越少。
# numericLabels = TRUE：表示模块将使用数值标签，返回数字而不是颜色作为模块的名字，后面可以再转换为颜色。
# pamRespectsDendro = FALSE：控制 PAM（Partitioning Around Medoids）聚类算法是否遵循树状结构。这里设置为 FALSE，表示不遵循树状结构。
# saveTOMs = F：用来控制是否保存共表达网络的拓扑重叠矩阵（TOM）。这里设置为 FALSE，表示不保存。最耗费时间的计算，有需要的话，大家存储起来，供后续使用。
# verbose = 3：用来控制输出的详细程度。设置为 3 可以输出较为详细的信息，用于跟踪函数的执行进程。
# corType：用来选择计算相关性的方法。默认“pearson”，还可以是“bicor”，“bicor”更能考虑离群点的影响bicor。

# 查看模块情况
table(net$colors) 
# 根据模块中基因数目的多少，降序排列，依次编号，比如 1 为最大模块，模块中基因最多。
# 0 表示没有分入任何模块的基因
# 使用层次聚类树展示各个模块，模块可视化

# 将基因模块的标签转换为对应的颜色编码
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# 绘制层次聚类树
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 这里用不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任何模块的那些基因，
# 如果灰色模块里面的基因太多，那么前期对表达矩阵挑选基因的步骤可能就不太合适。

########## 计算模块与表达量之间的关联性
# ME值，也就是获取eigengenes，每个ME代表一个模块
MEs <- net$MEs
head(MEs)[1:5, 1:5]
moduleLables <- net$colors
moduleColors <- labels2colors(net$colors)

geneTree <- net$dendrograms[[1]]

save(net, moduleLables, moduleColors, MEs, geneTree, file = "networkConstruction.RData")

# 将基因模块与性状进行关联

# 获取eigengenes，用颜色标签计算ME值
MEList <-  moduleEigengenes(exp_tpm, colors = moduleColors)
MEs0 <- MEList$eigengenes

# 查看用颜色标签计算的ME值
head(MEs0)[1:5, 1:5]

# 计算每个模块和每个性状之间的相关性
moduleTraitCor <- cor(MEs0, drug_screen_parp , use = "p");
head(moduleTraitCor)

# 计算显著性
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)

# 可视化相关性和P值
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3))

# 绘制热图
#pdf('labeledHeatmap.pdf', width = 8, height = 12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(drug_screen_parp),
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               # 如果前面的相关性较低，可以适当调整zlim值
               zlim = c(-1, 1), 
               main = paste("Module-trait relationships"))
dev.off()

# 在hypoxia性状中，MEdarkgrey正相关性最高，MEmagenta负相关性最高
# 温度性状中，MEgreenyellow正相关性最高，MEdarkred负相关性最高
# 获取模块名称
modNames <- substring(names(MEs0), 3)
modNames
# 计算模块与基因的相关性矩阵
geneModuleMembership <- as.data.frame(cor(exp_tpm, MEs0, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "");
geneModuleMembership[1:5, 1:5]

# 计算性状与基因的相关性矩阵 
# 只有连续型性状才能进行计算，如果是离散变量，在构建样本表时就转为0-1矩阵。
geneTraitSignificance <- as.data.frame(cor(exp_tpm, drug_screen_parp, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(drug_screen_parp), sep = "")
names(GSPvalue) <- paste("p.GS.", names(drug_screen_parp), sep = "")
head(geneTraitSignificance)


# 最后把两个相关性矩阵联合起来，指定感兴趣模块进行分析
module = "magenta"
pheno = "hypoxia"
modNames = substring(names(MEs0), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(drug_screen_parp))

# 获取模块内的基因
moduleGenes <- moduleColors == module

# 可视化基因与模块、表型的相关性，绘制散点图
par(mar = c(5, 5, 3, 3))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for LRG"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#### 提取指定模块的基因名
module = "greenyellow";
# Select module probes
probes = colnames(exp_tpm) ## 我们例子里面的probe就是基因
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
head(modProbes)
write.table(modProbes, paste0('Module_gene_name_', module, '_', pheno,'.txt'), row.names = F, quote = F)

# 如果使用WGCNA包自带的热图就很丑。
which.module="darkred";
dat=exp_tpm[,moduleColors==which.module ] 
plotMat(t(scale(dat)),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
exp_tpm[1:4,1:4]
dat=t(exp_tpm[,moduleColors==which.module ] )

library(pheatmap)
pheatmap(dat ,show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
group_list=drug_screen_parp$temperature
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =T,show_rownames = F)


########## 网络的可视化
# 计算非常漫长，建议只挑选部分基因
# 主要是可视化 TOM矩阵，WGCNA的标准配图
# 然后可视化不同 模块 的相关性 热图
# 不同模块的层次聚类图
# 还有模块诊断，主要是 intramodular connectivity
load('networkConstruction.RData')
load('WGCNA.RData')

if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  geneTree = net$dendrograms[[1]]; 
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6); 
  plotTOM = dissTOM^7; 
  diag(plotTOM) = NA; 
  #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 400
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  
  png("step7-Network-heatmap.png",width = 800,height = 600)
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  dev.off()
  
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  ## 只有连续型性状才能只有计算
  ## 这里把是否属 Luminal 表型这个变量0,1进行数值化
  Luminal = as.data.frame(design[,3]);
  names(Luminal) = "Luminal"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, Luminal))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  
  par(cex = 0.9)
  png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  dev.off()
  
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  ## 模块的进化树
  png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  dev.off()
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  ## 性状与模块热
  
  png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
  dev.off()
  
}
#这个非常消耗计算资源和时间，所以建议选取其中部分基因作图即可，我就没有画，而且根据下面的代码选取部分基因来作图！

#然后随机选取部分基因作图
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

png("step7-Network-heatmap.png",width = 800,height = 600)

TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
#最后画模块和性状的关系
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属 Luminal 表型这个变量0,1进行数值化
Luminal = as.data.frame(design[,3]);
names(Luminal) = "Luminal"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Luminal))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);

par(cex = 0.9)
png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
## 模块的进化树
png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热

png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


############ 模块的导出 ##########
# 主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape,VisANT等网络可视化软件。

# Recalculate topological overlap
library(WGCNA)
datExpr <- exp_tpm

TOM = TOMsimilarityFromExpr(datExpr, power = 12); 
save(TOM, file = 'TOM.Rdata')

load('TOM.Rdata')
load('networkConstruction.RData')
load('WGCNA.RData')


# Select module
module = "darkred";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

#首先是导出到VisANT

vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
#然后是导出到cytoscape

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);


#如果模块包含的基因太多，网络太复杂，还可以进行筛选，比如：

nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]


###########----------------CYTOSCAPE网络构建------------############
library(networkD3)
library(tidyverse)
library(vroom)
rm(list = ls())
gc()

module <- 'darkred'
df_node <- vroom::vroom(paste0('CytoscapeInput-nodes-',module,'.txt'))
#df_node <- cyt[[1]]
df_edge <- vroom::vroom(paste0('CytoscapeInput-edges-',module,'.txt'))
#df_edge <- cyt[[2]]

head(df_node)

### 数据过滤与筛选
df_edge <- df_edge %>% arrange(-weight) %>% head(100)
# 删除自身和自身相关位点
df_edge <- df_edge[which(df_edge$fromNode != df_edge$toNode),]
networkData <- df_edge[1:2]
p <- simpleNetwork(networkData,
                   linkDistance = 300,
                   opacity = 0.8,
                   fontSize = 10,
                   linkColour = '#dda8d8',
                   nodeColour = '#d608c2',
                   zoom = T)
p
saveNetwork(network = p,file = paste0('Network_', module, '.html'))

# 保存为pdf格式
library(webshot)
webshot('Network_yellowgreen.html','Network_yellowgreen.pdf')


### 高级绘图，可以添加组别、转录因子、表达量等
# 转换格式,把geneid转换为数字编号
df_edge_net <- df_edge[,c(1,2,3)] %>% as.data.frame()
df_node_net <- df_node[,c(1,3)] %>% as.data.frame()

colnames(df_edge_net) <- c("source" ,"target" ,"value")
colnames(df_node_net) <- c("name","group")

# 合并第一列和第二列，并取并集
merged_elements <- union_all(df_edge_net$source,df_edge_net$target) %>% unique()

# 对合并后的元素进行编号
element_numbers <- seq_along(merged_elements)

# 创建一个新的数据框，包含合并的元素和对应的编号
result_df <- data.frame(merged_elements, element_numbers)
result_df$element_numbers <- result_df$element_numbers-1


# 经过这一步处理后能够得到两个新的数据框，这就是绘制动态网络图的关键输入数据。在此基础上，我们还可以添加一些额外的信息，比如按照不同的分组将节点赋予不同的颜色，或者根据根据基因之间的正调控和负调控设置连接线的颜色。

# 生成模拟数据
df_edge_net$value <- c(runif(nrow(df_edge_net)/2,0,1),runif(nrow(df_edge_net)/2,0,5))
df_edge_net$color <- c(rep("red",50),rep("green",50))

#value值表示节点之间连线的权重大小，可以用来展示两个基因之间的关联程度，该值越大线越粗，关联性越强。

#color值可以用来设置连线的颜色，比如设置正调控为红色，负调控为绿色。

#除了设置节点与节点之间边的关系，还能设置单个节点的参数，比如通过下面的代码设置节点的大小用来表示基因的表达量，表达量高的基因节点直径越大。还可以用过Type将节点进行分组，比如转录因子为A组，目标基因为B组等等。

df_node_net <- result_df
df_node_net$size <- runif(nrow(df_node_net),0,20)
df_node_net$type <- rep(c("A","B","C"),10000)[1:nrow(df_node_net)]
colnames(df_node_net) <- c("name", "group", "size","type")



### 绘制动态网络图
# 使用映射表更新原始数据框的第一列和第二列
df_edge_net$source <- result_df$element_numbers[match(df_edge_net$source, result_df$merged_elements)]
df_edge_net$target <- result_df$element_numbers[match(df_edge_net$target, result_df$merged_elements)]

p <- forceNetwork(Links = df_edge_net, 
                  Nodes = df_node_net, 
                  Source = "source", 
                  Target = "target",
                  #linkColour=df_edge_net$color,
                  arrows=TRUE,
                  legend=TRUE,
                  Value = "value",
                  NodeID = "name",
                  Group = "group", 
                  bounded=F,
                  opacityNoHover = 0.5,
                  linkDistance = 100,
                  charge=-500,
                  #Nodesize='size',
                  # radiusCalculation = "Math.sqrt(d.nodesize,2)*5",
                  # linkWidth = JS("function(d) { return Math.sqrt(d.value)-4;}"),
                  # linkDistance=JS("function(d){return 1/(d.value)*100 }"),
                  opacity = 0.9,
                  zoom = T,
                  fontFamily = "Aril",
                  fontSize = 12) 
p
