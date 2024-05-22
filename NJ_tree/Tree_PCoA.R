# Get NJ tree and PCoA from the output of Calculate_pairwise_dist_simple.R
# Rscript this-code.R input_file output_file_prefix

args <- commandArgs(trailingOnly = TRUE)
input.file <- as.character(args[1])
output.prefix <- as.character(args[2])
output.file <- paste(output.prefix,"_tree_pcoa.rda",sep="")
output.pcoa <- paste(output.prefix,"_pcoa.pdf",sep="")
output.pcoa.txt <- paste(output.prefix,"_pcoa.txt",sep="")
output.tree.tre <- paste(output.prefix,"_tree.tre",sep="")
output.tree.nex <- paste(output.prefix,"_tree.nex",sep="")

load(file=input.file)
library(ape)

my.dist <- as.dist(distance.matrix)

my.pcoa <- pcoa(my.dist)
pdf(file=output.pcoa,width=6,height=6)
par(mar=c(5,5,1,1))
plot(my.pcoa$vectors[,1],my.pcoa$vectors[,2],pch=16,col=rgb(0,0,0,0.5),xlab="PCoA1",ylab="PCoA2")
dev.off()
write.table(my.pcoa$vectors,file=output.pcoa.txt,quote=F,col.names=T,row.names=F,sep="\t")

# Tree
attr(my.dist, "Labels") <- sample.name.all
my.tree <- nj(my.dist)
write.tree(my.tree, file=output.tree.tre)
write.nexus(my.tree, file=output.tree.nex)

save(sample.name.all,my.pcoa,my.tree,file=output.file)