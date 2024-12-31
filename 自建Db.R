rm(list = ls())
gc()

###制作可索引到物种的OrgDb包###
#下载和加载包#
BiocManager::install("AnnotationHub") 
BiocManager::install("AnnotationDbi")
BiocManager::install("rtracklayer")

library(AnnotationHub)
library(AnnotationDbi)
library(rtracklayer)

#索引与制作OrgDb#
hub <- AnnotationHub() #建立AnnotationHub对象保存到hub
#query(hub, 'Capra hircus')  #查询包含山羊(Capra hircus)的物种信息;结果有物种的各类信息需要进一步筛选
query(hub, 'largemouth bass')
#query(hub[hub$rdataclass == "OrgDb"] , "Capra hircus") #筛选我们需要OrgDb类型；也可将上一步与这一步合并成query(hub,'org.Capra hircus')进行搜索
#goat <- hub[['AH101444']]  #制作Capra hircus的OrgDb库；AH101444是Capra hircus对应的编号。
#goat　#查看goat
#help('select')


## 自建Db
egg <- fread('GO_database.txt',header = T)


## GO注释
gene_info <- egg %>%dplyr::select(GID = GeneID, GENENAME = GeneName) %>% na.omit() #把GID和GENENAME相应提取出来
goterms <- egg %>%dplyr::select(GeneID, Accession) %>% na.omit() %>% filter(str_detect(Accession,"GO")) #根据egg第一列和GO列的标题提取基因的GO注释。后面的`%>% filter(str_detect(GO,"GO"))`是筛选GO列值包含"GO"的行，删除空值。

# 根据逗号分隔GO号
library(stringr)
all_go_list <-  str_split(egg$Accession,",")
gene2go <- data.frame(GeneID = rep(egg$GeneID,
                                times = sapply(all_go_list, length)),
                      GO = unlist(all_go_list))
# 合并EVIDENCE号
gene_evidence <- egg[,c(1,4)]
gene2go <- merge(gene_info_2, gene_evidence)%>% na.omit()
colnames(gene2go) <- c('GID','EVIDENCE','GO')


#gene_info_2 <- unique(gene_info)  #我这里有行重复，你的没有的话可以忽略

#gene2go <- egg %>%dplyr::select(GID = GeneID, GO = Accession, EVIDENCE = EVIDENCE) %>% na.omit()  #同理把GID，GO和Evidence提取出来, Evidence随意什么都行。
#gene2go <- unique(gene2go)

write.table(gene_info_2,file="gene_info.txt",sep="\t",row.names=F,quote=F) 


#### 如果有KEGG注释的话运行以下代码

koterms <- egg %>%dplyr::select(GID = GeneID, KO=KO)%>%na.omit() #根据egg第一列和KEGG_ko列的标题提取基因的KEGG注释。

# 根据逗号分隔KO号
library(stringr)
all_ko_list <-  str_split(koterms$KO,",")
koterms_2 <- data.frame(GID = rep(koterms$GID,
                                       times = sapply(all_ko_list, length)),
                          KO = unlist(all_ko_list))
koterms_2 <- unique(koterms_2)


#### 对KEGG网站下载的json文件操作
if(!file.exists('kegg_info.RData')){
  
  library(jsonlite)
  library(purrr)
  library(RCurl)
  library(tibble)
  
  
    update_kegg <- function(json = "ko00001.json",file=NULL) {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    
    kegg <- fromJSON(json)
    
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
    
    save(pathway2name, ko2pathway, file = file)
  }
  
  update_kegg(json = "ko00001.json",file="kegg_info.RData")
  
}

load("kegg_info.RData")

# 在运行数据框合并前，需要做到两个数据框的列名是对应的，如有ko:需删除
library(stringr)
library(dplyr)
colnames(ko2pathway)=c("KO",'Pathway') #把ko2pathway的列名改为KO和Pathway，与koterms一致。
#koterms$KO <- str_replace_all(koterms_2$KO,"ko:","") #把koterms的KO值的前缀ko:去掉，与ko2pathway格式一致。
gene2pathway <- koterms_2 %>% left_join(ko2pathway, by = "KO") %>%dplyr::select(GID, Pathway) %>%na.omit() 
gene2pathway_uni <- unique(gene2pathway)
write.table(gene2pathway_uni,file="gene2pathway.txt",sep="\t",row.names=F,quote=F) 

##【optional】下面为可选步骤
gene2pathway_name<-left_join(gene2pathway,pathway2name,by="Pathway") #合并gene2pathway和pathway2name
write.table(gene2pathway_name,file="gene2pathway_name.txt",sep="\t",row.names=F,quote=F) #【optional】保存成文件gene2pathway_name.txt，可以给通用KEGG富集分析用。


#这里面的tax_id去NCBI的Taxonomy查Taxonomy ID信息：
#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
#BiocManager::install("AnnotationForge")
library(AnnotationForge)
gene2go <- gene2go[,c(1,3,2)]
makeOrgPackage(gene_info = gene_info, go = gene2go, ko = koterms_2,  pathway = gene2pathway_uni,
               version="0.0.1", maintainer='xiang<lixianghezi@163.com>',
               author='lixiang',outputDir=".", tax_id="489037", genus="Micropterus",
               species="Nigricans",goTable="go")


install.packages('org.MNigricans.eg.db',repos = NULL, type="source") #安装包

columns(org.MNigricans.eg.db)
head(keys(org.MNigricans.eg.db,keytype = 'GID'))
head(keys(org.MNigricans.eg.db,keytype = 'GO'))


## 测试数据
library(clusterProfiler)
library(data.table)
library(org.MNigricans.eg.db) #加载包

keytypes(org.MNigricans.eg.db)
gene <- fread("deg.txt",header = T)
gene_list <- gene$gene_name
#sig_gene <- gene[c(1:1000),1]
#gene_up <- gene[which(gene$change == c("UP")),1]  
#需要注意的是你的基因名和注释里的是否一样，如果不一样建议用SYMBOL更准确些
up_go <- enrichGO(gene_list,
                  OrgDb = org.MNigricans.eg.db,
                  keyType = 'GID',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  readable = F)
head(up_go)
up_go_bar <- barplot(up_go,showCategory=10,drop=T, title = "GO") #只显示前十个
up_go_bar




### KEGG富集
ko2gene <- fread('gene2pathway.txt')
ko2gene <- ko2gene[,c(2,1)]

ko2name <- fread('gene2pathway_name.txt')
ko2name <- ko2name[,c(2,3)]

KEGG_enrich <- enricher(gene_list,TERM2GENE = ko2gene ,TERM2NAME = ko2name,
                        pAdjustMethod = "BH",pvalueCutoff  = 0.05, qvalueCutoff  = 0.2)
KEGG_enrich_results <- KEGG_enrich@result

