gene.list <- read.csv("rlsset.sys.list.csv", header = T)

gene <- gene.list$Gene

pair <- c("Gene", "RLS")

for (i in 1:length(gene)) {
    gene.name <- as.character(gene[i])
    file.name <- paste("rlsset.renamed","/", gene[i], ".csv", sep="")
    rls.file <- read.csv(file.name, header = T)
    rls <- round(mean(rls.file$rls),5)
    pair <- rbind(pair, c(gene.name, rls))
}

write.table(pair, file="RLS.nonessential.mean.csv", sep=",", row.names=F, col.names=F)
