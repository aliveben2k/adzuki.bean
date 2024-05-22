
#check if packages are installed
#list.of.packages <- c("ggplot2", "grid", "gridExtra")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

library(ggplot2)

args <- commandArgs(trailingOnly = T)
filename <- args[1]
years_per_gen <- as.numeric(args[2])
colors <- c("dodgerblue","gold","red","navy", "orange3", "red4", "dodgerblue4", "tan1", "deeppink2", "purple1", "purple4")

#define x-axis ticks in ggplot
plot.ticks <- function(start, end) {
  tick.interval <- log10(2:9)
  tick.break <- c()
  start <- floor(start)
  end <- ceiling(log10(end))
  for (i in start:(end-1)){
    current.interval <- tick.interval + i
    tick.break <- c(tick.break, i, current.interval)
  }
  tick.break <- c(tick.break, end)
  return(as.numeric(tick.break))
}
#define x-axis marks in ggplot
plot.marks <- function(x) {
  y <- c()
  for (i in 1:length(x)){
    if(floor(x[i]) == x[i]){
      y <- c(y, parse(text=paste("10^",x[i])))
    } else {
      y <- c(y, "")
    }
  }
  return(y)
}

#read in population size
groups <- as.matrix(read.table(paste(filename, ".coal", sep = ""), nrow = 1))
time <- years_per_gen*t(as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = 1, nrow = 1)))
pop_size <- data.frame(time = numeric(0), pop_size = numeric(0), groups = numeric(0))
#pop_size <- data.frame(read.table(paste(filename, ".coal", sep = ""), skip = 2))
num_pops <- round(sqrt(dim(read.table(paste(filename, ".coal", sep = ""), skip = 2))[1]))
#num_pops <- (-1 + sqrt(1 + 8 * num_pops))/2

for(p1 in 1:num_pops){
  for(p2 in 1:p1){
    c        <- as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = (p1-1) * num_pops + p2 + 1, nrow = 1))[-c(1:2)]
    str      <- rep(paste(groups[p1]),length(c))
    if (groups[p1] == groups[p2]){
      pop_size <- rbind(pop_size, data.frame(time = time, pop_size = 0.5/c, groups = str))
    }
  }
}
#pop_size$time[which(pop_size$time > 1e5)] <- 1e5

#define number of colors in the plot
dis.cats <- length(groups)
if (dis.cats > 11){
  if(!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  colors <- getPalette(dis.cats)
} else {
  colors <- colors[1:dis.cats]
}

#plot
pop_size <- pop_size[which(pop_size$pop_size != "NaN"),]
tmax <- pop_size$time
xticks <- plot.ticks(0, max(tmax, na.rm = T))
yticks <- plot.ticks(0, max(pop_size$pop_size, na.rm = T))
options(scipen=999)

pop.plot <- ggplot(pop_size, aes(log10(time), log10(pop_size), color = groups)) +
  geom_step(linewidth = 0.5) +
  scale_x_continuous(breaks = xticks, labels = plot.marks(xticks), name = "Years ago") +
  scale_y_continuous(breaks = yticks, labels = plot.marks(yticks), name = expression(bolditalic(N[e]))) +
  scale_color_manual(aesthetics = "color", values = colors, name = "Genogroup") +
  theme_classic() + 
  theme(text = element_text(color = "black", size = 14),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)
  )

pdf(paste(filename,".pdf",sep = ""), width = 6, height = 3)
plot(pop.plot)
dev.off()


