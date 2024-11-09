library(plyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
# args[1]: matrix file
# args[2]: table file
# args[3]: path
# args[4]: private_SNP file

args <- commandArgs(trailingOnly = TRUE)
matrix <- read.table(args[1], sep = "\t")
pop.info <- read.table(args[2], header = F)
path <- args[3]
colnames(matrix) <- rownames(matrix)
rownames(pop.info) <- pop.info[,1]
mean_dist.pop.pi <- c()
y.max <- c()
pop <- unique(pop.info[,2])
for (i in 1:length(pop)){
  for (j in i:length(pop)){
    pop.sub <- pop.info[pop.info[,2] == pop[i],]
    pop.sub.2 <- pop.info[pop.info[,2] == pop[j],]
    data.sub <- matrix[pop.sub[,1], pop.sub.2[,1]]
    if (pop[i] == pop[j]){
      data.sub[lower.tri(data.sub)] <- NA
      data.sub[data.sub == 0] <- NA
    }
    data.sub <- as.vector(as.matrix(data.sub))
    data.mean <- mean(as.numeric(data.sub), na.rm = T)
    data.sd <- sd(as.numeric(data.sub), na.rm = T)
    if (pop[i] == pop[j]){
      pop.summary <- cbind(Genogroup = pop[i], m_dist = data.mean, m_distSD = data.sd)
    }
    mean_dist.pop.pi <- rbind(mean_dist.pop.pi, pop.summary)
  }
}
mean_dist.pop.pi <- unique(mean_dist.pop.pi)
privates <- read.table(args[4], sep = "\t", header = T)
private_flat <- melt(privates)

#points_private <- ddply(private_flat,.(variable),summarise,private_SNPs = mean(value), priSE = sqrt(var(value))/length(value))
#points_pi <- ddply(pi,.(GROUP),summarise,pi = mean(PI), piSE = sqrt(var(PI))/length(PI))
points_private <- ddply(private_flat,.(variable),summarise,private_SNPs = mean(value), priSD = sd(value))
colnames(points_private)[1] <- "Genogroup"
all_points <- merge(points_private, mean_dist.pop.pi, by = "Genogroup", all = TRUE)
colnames(all_points)[c(4,5)] <- c("m_dist", "m_distSD")
all_points[,c(4,5)] <- apply(all_points[,c(4,5)], 2, function(x) as.numeric(as.character(x)))

col.brewer.theme <- c("#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
pop <- unique(all_points$Genogroup)
col.pal <- col.brewer.theme[1:length(pop)]
all_points$Genogroup <- factor(all_points$Genogroup, levels=c("JP_LR","CN_LRN","CN_LRS","JP_WL","CN_WL","SOUTH_WL2","SOUTH_WL1"))

a <- ggplot(all_points, aes(x = m_dist,y = private_SNPs)) +
  geom_pointrange(data = all_points,aes(xmin = m_dist - m_distSD, xmax = m_dist + m_distSD,y = private_SNPs, colour = Genogroup), linewidth = 0.2, alpha = 0.5, size = 0) + 
  geom_pointrange(data = all_points,aes(ymin = private_SNPs - priSD, ymax = private_SNPs + priSD, x = m_dist,colour = Genogroup), linewidth = 0.2, alpha = 0.5, size = 0) +
  geom_point(aes(colour = Genogroup), size = 1, alpha = 1) +
  scale_color_manual(values=col.pal, name = "Group", breaks = c("JP_LR", "CN_LRN", "CN_LRS", "JP_WL", "CN_WL","SOUTH_WL2","SOUTH_WL1")) +
  theme_minimal() +
  theme(text = element_text(color = "black", size = 7.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title = element_text(color = "black", size = 7.5),
        #legend.position = "none",
        aspect.ratio = 1)
name_a <- "/private_mdist_plot.tiff"
tiff(paste0(path, name_a), units="cm", width = 8, height = 8, res = 600)
#pdf(paste0(path, "/output.pdf"),width=4,height=4)
print(a)
dev.off()
