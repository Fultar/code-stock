library(data.table)

rm(list = ls())
gc()

path <- c('../缺氧/36_results/')


### 提取显著性事件
for (splcing in c('RI','SE', 'MXE', 'A5SS', 'A3SS')){
  splicng_filename <- paste(path, splcing, '.MATS.JC.txt', sep = "")
  file_splicing <- fread(splicng_filename)
  file_splicing <- file_splicing[order(file_splicing$FDR),]
  file_splicing_FDR0.05 <- file_splicing[which(file_splicing$FDR < 0.05),]
  write.table(file_splicing_FDR0.05, paste0(path, splcing, '.MATS.JC.FDR0.05.txt'), 
              quote = F, row.names = F, sep = '\t')
}

### 合并事件

for (path  in c('../缺氧/36_results/','../缺氧/32_results/','../缺氧/28_results/',
                '../温度/28VS32/','../温度/28VS36/','../温度/32VS36/')){
  file_A3SS <- fread(paste0(path, 'A3SS.MATS.JC.FDR0.05.txt'))
  file_A3SS$AS_events <- paste('A3SS_',file_A3SS$GeneID,'_',file_A3SS$chr, '_', file_A3SS$shortES,'_',file_A3SS$shortEE)
  
  file_A5SS <- fread(paste0(path, 'A5SS.MATS.JC.FDR0.05.txt'))
  file_A5SS$AS_events <- paste('A5SS_',file_A5SS$GeneID,'_',file_A5SS$chr, '_', file_A5SS$shortES,'_',file_A5SS$shortEE)
  
  file_RI <- fread(paste0(path, 'RI.MATS.JC.FDR0.05.txt'))
  file_RI$AS_events <- paste('RI_',file_RI$GeneID,'_',file_RI$chr, '_', file_RI$riExonStart_0base,'_',file_RI$riExonEnd)
  
  file_MXE <- fread(paste0(path, 'MXE.MATS.JC.FDR0.05.txt'))
  file_MXE$AS_events <- paste('MXE_',file_MXE$GeneID,'_1st_',file_MXE$chr, '_', file_MXE$'1stExonStart_0base','_',file_MXE$'1stExonEnd',
                              '_2nd_',file_MXE$'2ndExonStart_0base','_',file_MXE$'2ndExonEnd')
  
  file_SE <- fread(paste0(path, 'SE.MATS.JC.FDR0.05.txt'))
  file_SE$AS_events <- paste('SE_',file_SE$GeneID, '_', file_SE$chr, '_', file_SE$exonStart_0base,'_',file_SE$exonEnd)
  
  file_merge <- rbind(file_A3SS[,c(-6,-7,-8,-9,-10,-11)],
                      file_A5SS[,c(-6,-7,-8,-9,-10,-11)],
                      file_MXE[,c(-6,-7,-8,-9,-10,-11,-12,-13)],
                      file_RI[,c(-6,-7,-8,-9,-10,-11)],
                      file_SE[,c(-6,-7,-8,-9,-10,-11)])
  
  write.table(file_merge, file = paste0('DAS_', strsplit(path, '/')[[1]][3], '.txt'),
              sep = '\t', quote = F, row.names = F)
}



### 旧代码
for (path  in c('../缺氧/36_results/','../缺氧/32_results/','../缺氧/28_results/',
                '../温度/28VS32/','../温度/28VS36/','../温度/32VS36/')){
  file_A3SS <- fread(paste0(path, 'A3SS.MATS.JC.FDR0.05.txt'))
  file_A3SS$AS_events <- paste('A3SS_',file_A3SS$GeneID,'_',file_A3SS$shortES,'_',file_A3SS$shortEE)
  
  file_A5SS <- fread(paste0(path, 'A5SS.MATS.JC.FDR0.05.txt'))
  file_A5SS$AS_events <- paste('A5SS_',file_A5SS$GeneID,'_',file_A5SS$shortES,'_',file_A5SS$shortEE)
  
  file_RI <- fread(paste0(path, 'RI.MATS.JC.FDR0.05.txt'))
  file_RI$AS_events <- paste('RI_',file_RI$GeneID,'_',file_RI$riExonStart_0base,'_',file_RI$riExonEnd)
  
  file_MXE <- fread(paste0(path, 'MXE.MATS.JC.FDR0.05.txt'))
  file_MXE$AS_events <- paste('MXE_',file_MXE$GeneID,'_1st_',file_MXE$'1stExonStart_0base','_',file_MXE$'1stExonEnd',
                              '_2nd_',file_MXE$'2ndExonStart_0base','_',file_MXE$'2ndExonEnd')
  
  file_SE <- fread(paste0(path, 'SE.MATS.JC.FDR0.05.txt'))
  file_SE$AS_events <- paste('SE_',file_SE$exonStart_0base,'_',file_SE$exonEnd)
  
  file_merge <- rbind(file_A3SS[,c(-6,-7,-8,-9,-10,-11)],
                      file_A5SS[,c(-6,-7,-8,-9,-10,-11)],
                      file_MXE[,c(-6,-7,-8,-9,-10,-11,-12,-13)],
                      file_RI[,c(-6,-7,-8,-9,-10,-11)],
                      file_SE[,c(-6,-7,-8,-9,-10,-11)])
  
  write.table(file_merge,paste0(path, 'ALL_DAS.MATS.JC.txt'), sep = '\t', quote = F, row.names = F)
}
