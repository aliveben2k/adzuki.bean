# args[1]: matrix file
# args[2]: table file
# args[3]: path
library(tidyverse)
library(ggplot2)
library(ape)
library(reshape2)
#library(RColorBrewer)
#library(viridisLite)
#col.brewer.theme <- "Paired" #for RColorBrewer
#col.brewer.theme <- "H" #for viridisLite

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

#col.brewer.theme <- c("dodgerblue","orange","red","navy", "red4", "gold", "darkviolet", "maroon1", "peru", "coral1") #user_defined
#col.brewer.theme <- c("black","grey50","grey75","dodgerblue", "orange", "red", "darkviolet", "maroon1", "peru", "coral1") #user_defined
col.brewer.theme <- c("black","grey50","grey75","#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")

args <- commandArgs(trailingOnly = TRUE)
matrix <- read.table(args[1], sep = "\t")
group.info <- read.table(args[2], header = T)
colnames(matrix) <- rownames(matrix)
rownames(group.info) <- group.info[,1]
matrix <- merge(group.info, matrix, by = 0, all.x = F)
group.info <- group.info[c(matrix[,1]),]
matrix <- matrix[,2:ncol(matrix)]
colnames(matrix)[2] <- "group"
matrix <- matrix[which(matrix[,1] != "DRR199645"),which(colnames(matrix) != "DRR199645")]
matrix <- matrix[, 3:ncol(matrix)]
rownames(matrix) <- colnames(matrix)

uni.group <- unique(group.info[,2])
uni.group <- uni.group[grepl("JP_WL|JP_LR|CN_LRN", uni.group, ignore.case = T)]

final.table <- c()
y.max <- c()
for(i in 1:length(uni.group)){
  for (k in i:length(uni.group)){
    ##original
    #p1.subset <- matrix[which(matrix$group == uni.group[i]),]
    #p1.subset.mean <- colMeans(p1.subset[,3:ncol(p1.subset)])
    #message(uni.group[k])
    #p2.idx <- grep(uni.group[k], matrix$group)
    #p1.vs.p2.dist <- p1.subset.mean[p2.idx]
    #message(paste0(uni.group[i], " vs. ",uni.group[k]))
    ##CRLee
    p1.vs.p2.dist <- as.matrix(matrix[which(group.info[,2] == uni.group[i]),which(group.info[,2] == uni.group[k])])
    d <- density(as.numeric(p1.vs.p2.dist[upper.tri(p1.vs.p2.dist)]))
    d <- data.frame(x = d$x, y = d$y)
    y.max <- c(y.max, max(d$y))
    if (uni.group[i] == uni.group[k]){
      sub.group.name <- rep(uni.group[i], nrow(d))
    } else {
      sub.group.name <- rep(paste0(uni.group[i], " vs. ", uni.group[k]), nrow(d))
    }
    d <- cbind(group = sub.group.name, d)
    final.table <- rbind(final.table, d)
  }
}
uni.group <- unique(final.table[,1])
col.pal <- col.brewer.theme[1:length(uni.group)] #user_defined

x.max <- max(as.numeric(final.table[,2])) + max(as.numeric(final.table[,2])) * 0.15
y.max <- max(y.max)

col.brewer.theme <- c("#255271","#4b8bcb","gold","#f9bbb9", "#c6b1d4", "#90a7b7")
final.table$group <- factor(final.table$group, levels = c("JP_WL","JP_LR","CN_LRN","JP_LR vs. JP_WL","JP_WL vs. CN_LRN","JP_LR vs. CN_LRN"))

out.plot <- ggplot(data.frame(final.table), aes(x = as.numeric(x), y = as.numeric(y), group = group)) +
  #stat_density(aes(x=as.numeric(x)), geom="line", position="identity", size = 1) +
  #geom_line(position="identity", size = 1) +
  geom_area(position="identity", aes(fill = group), size = NULL, alpha = 0.8) +
  scale_color_manual(values = col.pal) +
  scale_fill_manual(values = col.pal, breaks = c("JP_WL","JP_LR","CN_LRN","JP_LR vs. JP_WL","JP_WL vs. CN_LRN","JP_LR vs. CN_LRN")) +
  scale_x_continuous(breaks = seq(0.000, 0.025, by = 0.01), name = "Pairwise distance") +
  ylab("Density") +
  theme_minimal() +
  theme(text = element_text(color = "black", size = 7.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title = element_text(color = "black", size = 7.5),
        legend.position = "none",
        aspect.ratio = 1)
  #xlim(0, x.max)# + ylim(0,y.max*1.3)
  

pdf(paste0(args[3],"/distance_density.pdf"))
print(out.plot)
dev.off()

name_a <- "E:\\cloud_files\\Lee\ lab\\Vigna\ angularis\\new_gatk4\\final_no_out\\20211222_0.7\\distance_density_north.tiff"
tiff(name_a, units="cm", width = 4, height = 4, res = 600)
print(out.plot)
dev.off()
