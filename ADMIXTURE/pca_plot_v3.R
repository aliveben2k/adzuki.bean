# args[1]: matrix file
# args[2]: table file
# args[3]: path
# args[4]: admix group true or false
library(tidyverse)
library(ggplot2)
library(ape)
#library(RColorBrewer)
#library(viridisLite)
#col.brewer.theme <- "Paired" #for RColorBrewer
#col.brewer.theme <- "H" #for viridisLite
col.brewer.theme <- c("dodgerblue","orange","red","navy", "red4", "gold", "darkviolet", "maroon1", "peru", "coral1") #user_defined
args <- commandArgs(trailingOnly = TRUE)
matrix <- read.table(args[1], sep = "\t")
my.dist <- as.dist(matrix)
my.pcoa <- pcoa(my.dist)
eigenval <- my.pcoa[["values"]][["Eigenvalues"]][1:20]
eigenvec <- as.data.frame(my.pcoa[["vectors"]][,1:20])
names(eigenvec)[1:ncol(eigenvec)] <- paste0("PC", 1:ncol(eigenvec))

pca <- read.table(args[2])
rownames(pca) <- pca[,1]
pca <- pca[,2:3]
#names(pca)[1] <- as.character("id")
names(pca)[1] <- as.character("group")
names(pca)[2] <- as.character("type")
pca <- merge(pca, eigenvec, by = 0, all.x = F)

types.vec <- unique(pca$type)
type <- rep(NA, length(pca$type))
for (y in 1:length(types.vec)) {
  type[grep(types.vec[y], pca$type)] <- types.vec[y]
}
group.vec <- unique(pca$group)
group <- rep(NA, length(pca$group))
for (x in 1:length(group.vec)) {
  group[grep(group.vec[x], pca$group)] <- group.vec[x]
}
group_type <- paste0(group, "_", type)
pca <- as_tibble(data.frame(pca, group, type, group_type))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
colnames(pca)[2] <- "Groups"
colnames(pca)[3] <- "Types"
group.break <- c()
if (as.numeric(args[4]) == 1){
  admix <- as.numeric(grep("admix", group.vec))
  for (y in 1:length(group.vec)){
    if (y != admix){
      group.break <- c(group.break, group.vec[y])
    }
  }
  group.break <- sort(group.break)
  group.break <- c(group.break, group.vec[admix])
#  col.pal <- brewer.pal(n = length(group.vec)-1, name = col.brewer.theme) #for RColorBrewer
#  col.pal <- viridis(length(group.vec)-1, alpha = 1, begin = 0, end = 1, option = col.brewer.theme) #for viridisLite
  col.pal <- col.brewer.theme[1:length(group.vec)-1] #user_defined
  col.pal <- c(col.pal, "grey70")
} else {
  group.break <- sort(group.vec)
#  col.pal <- brewer.pal(n = length(group.vec), name = col.brewer.theme) #for RColorBrewer#
#  col.pal <- viridis(length(group.vec), alpha = 1, begin = 0, end = 1, option = col.brewer.theme) #for viridisLite
  col.pal <- col.brewer.theme[1:length(group.vec)] #user_defined
}
if (length(na.omit(pca$Types)) == 0){
  b <- ggplot(pca, aes(PC1, PC2, color = Groups)) + geom_point(alpha = 0.8, size = 1.5) +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = col.pal, breaks = group.break) +
    theme_minimal() + theme(aspect.ratio = 1, 
                            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                            panel.grid = element_blank(),
                            axis.ticks = element_line(colour = "black", size = 0.5),
                            axis.title = element_text(face = "bold"),
                            axis.text = element_text(color = "black", size = 14),
                            #axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.title = element_text(face = "bold"),
                            text = element_text(color = "black", size = 16)) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
} else {
  b <- ggplot(pca, aes(PC1, PC2, color = Groups, shape = Types)) + geom_point(alpha = 0.8, size = 1.5) +
    scale_fill_manual(aesthetics = c("colour", "fill"), values = col.pal, breaks = group.break) +
    theme_minimal() + theme(aspect.ratio = 1, 
                            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
                            panel.grid = element_blank(),
                            axis.ticks = element_line(colour = "black", size = 0.5),
                            axis.title = element_text(face = "bold"),
                            axis.text = element_text(color = "black", size = 14),
                            #axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.title = element_text(face = "bold"),
                            text = element_text(color = "black", size = 16)) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
}
if (!is.na(args[3])){
  name_a <- paste0(args[3],"/PCA_analysis.pdf")
}
name_b <- paste0(args[3],"/PCA_K",length(group.vec)-as.numeric(args[4]),".pdf")
write.table(pca, file = paste0(args[3], "/PCA_K",length(group.vec)-as.numeric(args[4]), "_table.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
save(pve, pca, file = paste0(args[3], "/PCA_K",length(group.vec)-as.numeric(args[4]), "_table.rda"))
pdf(name_a)
print(a)
dev.off()
pdf(name_b, width = 5.5, height = 5.5)
print(b)
dev.off()
