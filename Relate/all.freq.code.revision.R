library(Rcpp)
library(RcppCNPy)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
path <- c(); color.series <- c()
for (i in 1:length(args)){
  if (args[i] == '-p'){
    path <- as.character(args[i+1])
    if (grepl("/$|\\\\$", path)){
      path <- sub("/$|\\\\$",'',path)
    }
  }
  if (args[i] == '-c'){
    color.series <- as.character(args[i+1])
  }
}

if (length(color.series) == 0){
  color.series <- "YlGnBu"
}

files <- c(); seris <- c()
file.loci <- c()
files <- list.files(path, pattern="npy$", full.names=TRUE)
for (i in 1:length(files)){
  path.tmp <- unlist(strsplit(files[i],"/|\\\\"))
  seris <- path.tmp[length(path.tmp)]
  file.locus <- unlist(strsplit(files[i], "_CLUES_"))[2]
  file.locus <- sub("\\.epochs\\.npy$|\\.freqs\\.npy$|\\.post\\.npy$","", file.locus) #no path
  file.loci <- c(file.loci, file.locus)
}
file.loci <- unique(file.loci)

plot.files <- c()
y.max.broken <- list()
for (i in 1:length(file.loci)){
  files.curr <- files[grepl(file.loci[i], files)]
  file.prefix <- gsub("\\.epochs\\.npy$|\\.freqs\\.npy$|\\.post\\.npy$","", files.curr) #with path
  file.prefix <- unique(file.prefix)
  post.exp.locus.all <- c()
  freqs.locus.all <- c()
  y.max.broken.curr <- c()
  for (j in 1:length(file.prefix)){
    epochs <- npyLoad(paste0(file.prefix[j],".epochs.npy"))
    freqs <- npyLoad(paste0(file.prefix[j], ".freqs.npy"))
    post <- npyLoad(paste0(file.prefix[j], ".post.npy"))
    colnames(post) <- epochs[1:length(epochs)-1]
    post.exp <- exp(post)
    if (j == 1){
      post.exp.locus.all <- post.exp
      freqs.locus.all <- freqs 
    } else {
      post.exp.locus.all <- post.exp.locus.all + post.exp
      freqs.locus.all <- freqs.locus.all + freqs 
    }
    y.max <- apply(post.exp, 2, which.max)
    y.max <- y.max/150
    y.max <- as.data.frame(cbind(epochs[1:length(epochs)-1], y.max))
    colnames(y.max) <- c("x","y")
    if (j == 1){
      locus.name <- rep(file.loci[i], nrow(y.max))
      y.max.broken.curr <- data.frame(y.max)
      y.max.broken.curr <- cbind(locus.name, y.max.broken.curr)
      colnames(y.max.broken.curr) <- c("locus", "time", j)
    } else {
      y.max.broken.curr <- cbind(y.max.broken.curr, y.max$y)
      colnames(y.max.broken.curr)[ncol(y.max.broken.curr)] <- j
    }
  }
  y.max.broken[[i]] <- y.max.broken.curr
  post.exp.locus.all <- post.exp.locus.all / length(file.prefix)
  freqs.locus.all <- freqs.locus.all / length(file.prefix)
  time <- epochs[1:length(epochs)-1]
  if (!file.exists(paste0(path, "/For_plot_", file.loci[i], ".rda"))){
    save(post.exp.locus.all, freqs.locus.all, time, color.series, file = paste0(path, "/For_plot_", file.loci[i], ".rda"))
  }
  plot.files <- c(plot.files, paste0(path, "/For_plot_", file.loci[i], ".rda"))
}

#plotting all loci in one figure
y.max.all <- c()
if (!file.exists(paste0(path, "/plot_data_all_loci.rda"))){
  for (i in 1:length(plot.files)){
    load(plot.files[i])
    y.max <- apply(post.exp.locus.all, 2, which.max)
    y.max <- y.max/150
    y.max <- as.data.frame(cbind(time, y.max))
    colnames(y.max) <- c("x","y")
    post.exp.locus.all <- melt(post.exp.locus.all)
    post.exp.locus.all$Var1 <- post.exp.locus.all$Var1/max(post.exp.locus.all$Var1)
    colnames(post.exp.locus.all) <- c("Frequency", "Time", "value")
    label.curr <- rep(file.loci[i], nrow(y.max))
    y.max <- cbind(y.max, label.curr)
    colnames(y.max)[3] <- "label"
    y.max$freq.bottom <- NA
    y.max$freq.top <- NA
    for (i in 1:nrow(y.max)) {
      time <- y.max$x[i]
      freq.max.prob <- y.max$y[i]
      post.exp.this.time <- post.exp.locus.all[post.exp.locus.all$Time == time,]
      # Which of the 14999 row has max post.prob
      which.row <- which(post.exp.this.time$Frequency == freq.max.prob)
      # From row 1 to this row -1, -2, which row has sum of prob >= 0.025
      which.row.bottom <- c()
      for (how.many.rows.back in 1:which.row) {
        if (sum(post.exp.this.time$value[1:(which.row - how.many.rows.back)]) <= 0.025) {
          which.row.bottom <- which.row - how.many.rows.back
          break
        }
      }
      # From this row +1, +2, to the last row, which row has sum of prob >= 0.025
      which.row.top <- c()
      last.row <- nrow(post.exp.this.time)
      for (how.many.rows.forward in 1:(last.row-which.row)) {
        if (sum(post.exp.this.time$value[last.row:(which.row + how.many.rows.forward)]) <= 0.025) {
          which.row.top <- which.row + how.many.rows.forward
          break
        }
      }
      # Get the frequency
      freq.bottom <- post.exp.this.time$Frequency[which.row.bottom]
      freq.top <- post.exp.this.time$Frequency[which.row.top]
      # Fill in the table y.max
      y.max$freq.bottom[i] <- freq.bottom
      y.max$freq.top[i] <- freq.top
    }
    y.max.all <- rbind(y.max.all, y.max)
  }
  colnames(y.max.all) <- c("Time","Frequency","label","freq.bottom","freq.top")
  save(y.max.all, file = paste0(path, "/plot_data_all_loci.rda"))
}

#adzuki bean use
plot.data <- data.frame(label = character(), Time = numeric(), Frequency = numeric(), freq.uppder = numeric, freq.lower = numeric())
for (i in 1:length(y.max.broken)){ #loci list
  y.max.curr <- y.max.broken[[i]]
  name.locus <- unique(y.max.curr[,1])
  #three.lines <- c()
  plot.data.tmp <- c()
  for (j in 1:nrow(y.max.curr)){
    point.median <- median(as.numeric(y.max.curr[j,3:ncol(y.max.curr)])) #median
    point.upper <- sort(as.numeric(y.max.curr[j,3:ncol(y.max.curr)]))[ncol(y.max.curr)-2-3] # the fourth largest
    point.lower <- sort(as.numeric(y.max.curr[j,3:ncol(y.max.curr)]))[1+3] # the fourth smallest
    point.values <- c(name.locus, y.max.curr[j,2], point.median, point.upper, point.lower)
    plot.data <- rbind(plot.data, point.values)
  }
}
colnames(plot.data) <- c("label","Time","Frequency","freq.upper","freq.lower")
plot.data[,2:ncol(plot.data)] <- plot.data[,2:ncol(plot.data)] %>% mutate_if(is.character, as.numeric)
save(y.max.broken, plot.data, file = paste0(path, "/final_all_loci.broken.stick.rda"))
#end of adzuki bean use
plot.data$label <- factor(plot.data$label, levels = c("CN_LRNCN_LRSJP_LRSOUTH_LR_chr1_13815643-13815643", "CN_LRNCN_LRSJP_LRSOUTH_LR_chr4_4148520-4148520", "CN_LRNCN_LRSJP_LRSOUTH_LR_chr7_2492619-2492619"))
comp <- ggplot() +
  geom_ribbon(data = plot.data, aes(x = Time, ymin = freq.lower, ymax = freq.upper, fill = label), alpha = 0.2) +
  geom_line(data = plot.data, aes(x = Time, y = Frequency, color = label), linewidth = 0.5) +
#adzuki bean use
  geom_segment(mapping= aes(x=16780.2, xend=20000, y=0.95, yend=0.95), linewidth = 0.5, color = "deepskyblue3", inherit.aes = F) + #ANR1
  geom_segment(mapping= aes(x=10740.8, xend=13260.8, y=0.9, yend=0.9), linewidth = 0.5, color = "darkgoldenrod3", inherit.aes = F) + #PAP1
  geom_segment(mapping= aes(x=7512.51, xend=10289.9, y=0.85, yend=0.85), linewidth = 0.5, color = "darkorchid3", inherit.aes = F) + #MYB26
#end of adzuki bean use
  scale_color_manual(values = c("deepskyblue3","darkgoldenrod3","darkorchid3"),
                     name = "Gene", labels = c("ANR1","PAP","MYB26")) +
  scale_fill_manual(values = c(c("deepskyblue3","darkgoldenrod3","darkorchid3")),
                    name = "95% CI", labels = c("ANR1","PAP","MYB26")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,20000), breaks = seq(0, 18000, 9000)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        #axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        text = element_text(size = 7.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(-0.2, 'cm'),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title.x = element_text(color = "black", size = 7.5),
        axis.title.y = element_text(color = "black", size = 7.5),
        legend.key.size = unit(0.3, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank()
  )
tiff_out = paste0(path, "/final_all_loci.median.tiff")
tiff(tiff_out, units = "cm", res = 600, width = 8, height = 4)
print(comp)
dev.off()



