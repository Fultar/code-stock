library(data.table)


file <- fread('Integrated_Function.annotation.txt', header = T)

string

result <- regexpr(pattern= ".GO:{7}.", text = file$GO_annotation)

start <- attr(result,"capture.start")
length <- attr(result,"capture.length")
name <- attr(result,"capture.name")
geneID <- ifelse(start > 0, 
                 substr(file$GO_annotation, start[,name],start[,name] + length[,name]-1),NA)
