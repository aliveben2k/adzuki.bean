library(Rcpp)
library(RcppCNPy)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
file.name <- c(); color.series <- c()
for (i in 1:length(args)){
  if (args[i] == '-f'){
    file.name <- as.character(args[i+1])
  }
  if (args[i] == '-c'){
    color.series <- as.character(args[i+1])
  }
}

if (length(color.series) == 0){
  color.series <- "YlGnBu"
}

epochs <- npyLoad(paste0(file.name, ".epochs.npy"))
freqs <- npyLoad(paste0(file.name, ".freqs.npy"))
post <- npyLoad(paste0(file.name, ".post.npy"))
out <- paste0(file.name, ".pdf")
rda <- paste0(file.name, ".rda")
colnames(post) <- epochs[1:length(epochs)-1]
post.exp <- exp(post)
y.max <- c()
for (i in 1:length(colnames(post.exp))){
  y.max.curr.val <- 0
  for (j in 1:length(post.exp[,1])){
    if (post.exp[j,i] > y.max.curr.val){
      y.max.curr.val <- post.exp[j,i]
      y.max.curr <- j
    }
  }
  y.max <- c(y.max, y.max.curr)
}
y.max <- y.max/150
y.max <- as.data.frame(cbind(epochs[1:length(epochs)-1], y.max))
colnames(y.max) <- c("x","y")
post.exp <- melt(post.exp)
post.exp$Var1 <- post.exp$Var1/max(post.exp$Var1)
colnames(post.exp) <- c("Frequency", "Time", "value")
fig_width <- 5
fig_height <- 3
save(post.exp, y.max, color.series, fig_width, fig_height, file = rda)

fig <- ggplot(post.exp, aes(x=Time, y=Frequency, fill=value)) + 
  geom_tile() +
  geom_line(data = y.max, aes(x = x, y = as.numeric(y)), color = "black", inherit.aes = FALSE, linewidth = 0.2) +
  #geom_smooth(data = y.max, aes(x = x, y = as.numeric(y)), color = "black", method = "gam", inherit.aes = FALSE, linewidth = 0.5) +
  scale_fill_distiller(palette = color.series, 
                       direction = 1,
                       name = "Posterior probability",
                       guide = guide_colourbar(title.position = "right", 
                                               barheight = unit(fig_height*0.8, "inches"), 
                                               barwidth = unit(fig_width*0.01, "inches"),
                                               ticks = F,
                                               ticks.color = "black",
                                               ticks.linewidth = 1,
                                               label.vjust = 0.5, 
                                               title.hjust = 0.5)
                       ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        text = element_text(face = "bold", size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0)),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0), angle = -90),
        legend.text = element_text(size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0)),
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5))

pdf(out, width=fig_width, height=fig_height)
print(fig)
dev.off()