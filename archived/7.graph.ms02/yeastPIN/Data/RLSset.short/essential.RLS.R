gene.list <- read.csv("rlsset.sys.list.csv", header = T)
gene <- gene.list$Gene

essential <- read.csv("essential.list", header=T)
ess.gene <- essential$Gene

pair <- c("Gene", "RLS")

for (i in 1:length(ess.gene)) {
    gene.name <- as.character(ess.gene[i])
    test <- which(gene %in% gene.name)
    if (length(test) == 0) {
       RLS <- 0.0
    } else {
       RLS <- NA 
    }
    pair <- rbind(pair, c(gene.name, RLS))
}

write.table(pair, file="RLS.essential.mean.csv", sep=",", row.names=F, col.names=F)
