options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

filename <- gsub(" ", "", args[1])
rm(args)

library("cluster")
library("ape")

dat <- read.csv(file=filename, head=F)
dv <- diana(dat, diss=T)
num <- ncol(dat)
clusters <- cutree(as.hclust(dv), k=1:num)
tree <- as.phylo(as.hclust(dv))

newick_name <- paste(filename, ".newick", sep = "")
write.tree(phy=tree, file=newick_name)
clusters_name <- paste(filename, ".clusters.csv", sep = "")
write.csv(clusters, file=clusters_name)

