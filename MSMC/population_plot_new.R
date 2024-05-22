if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
#library(ggalt)
os <- Sys.info()[['sysname']]
colors <- c("#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")

args <- commandArgs(trailingOnly = TRUE)
mu <- as.numeric(args[1]) #mutation rate
mu_out <- as.character(mu)
gen <- as.numeric(args[2]) #generation time
random <- as.character(args[3]) #serial number/identifier
path <- sub("/$", "", args[4]) #folder path containing input files

dirs <- list.dirs(path, full.names = T, recursive = T)
files <- c()
for (k in 1:length(dirs)){
  files.tmp <- list.files(dirs[k], pattern="all\\.final\\.out$", full.names=TRUE)
  files.tmp <- files.tmp[!grepl("combined", files.tmp)]
  files <- c(files, files.tmp)
}
if (length(files) == 0){
  for (k in 1:length(dirs)){
    files.tmp <- list.files(dirs[k], pattern="final\\.txt$", full.names=TRUE)
    files.tmp <- files.tmp[!grepl("combined", files.tmp)]
    tmp <- unlist(strsplit(files.tmp[1],split="/"))
    folder.name <- tmp[length(tmp)-1]
    files <- c(files, files.tmp[!grepl(paste0(folder.name, "_"), files.tmp)])
  }
}
if (length(files) == 0){
  message("Cannot find the file for plotting.")
  quit("no")
}

#combine all population table together
final.table <- c()
for (i in 1:length(files)){
  file_names <- unlist(strsplit(files[i],"/", fixed=T))
  pop_name <- unlist(strsplit(file_names[length(file_names)],"_|\\."))
  pop_name <- paste0(pop_name[1],"_",pop_name[2])
  possibleError <- tryCatch({
  curr.pop.table <- read.table(files[i], header=T, sep="\t")}, 
  error = function(e) e)
  if(inherits(possibleError, "error")){ 
    message(paste0(files[i], ' is empty. Skip.'))
    next
  }
  pop <- rep(pop_name, nrow(curr.pop.table))
  x.left <- log10(curr.pop.table$mean.left/mu*gen)
  x.center <- (log10(curr.pop.table$mean.left/mu*gen) + log10(curr.pop.table$mean.right/mu*gen))/2
  x.right <- log10(curr.pop.table$mean.right/mu*gen)
  y <- log10((1/curr.pop.table$mean.lambda)/(2*mu))
  ci.upper.log10 <- log10((1/curr.pop.table$upper.bound)/(2*mu))
  ci.lower.log10 <- log10((1/curr.pop.table$lower.bound)/(2*mu))
  curr.pop.table <- cbind(curr.pop.table, x.left, x.center, x.right, y, ci.upper.log10, ci.lower.log10, pop)
  final.table <- rbind(final.table, curr.pop.table)
}

#define number of colors in the plot
final.table[,ncol(final.table)] <- as.character(final.table[,ncol(final.table)])
dis.cats <- length(unique(final.table[,ncol(final.table)]))
if (dis.cats > 11){
  if(!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  colors <- getPalette(dis.cats)
} else {
  colors <- colors[1:dis.cats]
}

#define x-axis ticks in ggplot
plot.ticks <- function(start, end) {
  tick.interval <- log10(2:9)
  tick.break <- c()
  start <- floor(start)
  end <- ceiling(end)
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

##optional code
##define the right boundary of the x-axis for each population
#V. angularis use
# all.names <- c("CN_LRN","CN_WL","JP_LR","JP_WL","SOUTH_WL1","SOUTH_WL2")
# all.limit <- c(50000, 140000, 50000, 80000, 1000000, 1000000)
# all.limit <- log10(all.limit)
# subsets <- c()
# for (i in 1:length(all.names)){
#   subset.tmp <- final.table[final.table$pop == all.names[i],]
#   subset.tmp <- subset.tmp[subset.tmp$x.left <= all.limit[i],]
#   subsets <- rbind(subsets, subset.tmp)
# }
# if (length(subsets) > 0){
#   final.table <- subsets
# }
##define population order on the figure
#V. angularis use
#colors <- c("dodgerblue","gold","navy", "orange3", "red4", "dodgerblue4", "tan1", "deeppink2", "purple1", "purple4")
#colors <- c("#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
#colors <- c("#4b8bcb","gold","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
#final.table$pop <- factor(final.table$pop, levels = c("JP_LR", "CN_LRN", "JP_WL", "CN_WL", "SOUTH_WL2","SOUTH_WL1"))

#ggplot
tmax <- final.table$x.left[which(final.table$x.left < Inf)]
xticks <- plot.ticks(0, max(tmax, na.rm = T))
yticks <- plot.ticks(0, max(final.table$y, na.rm = T))
options(scipen=999)
pop.plot <- ggplot(final.table) +
  #geom_ribbon(aes(x = x.left, ymin = ci.lower.log10, ymax = ci.upper.log10, fill = pop), alpha = 0.5) +
  geom_rect(aes(xmin = x.left, xmax = x.right, ymin = ci.lower.log10, ymax = ci.upper.log10, fill = pop), alpha = 0.5) +
  geom_step(linewidth = 0.2, aes(x.left, y, color = pop)) +
  #geom_line(aes(x.left, y, color = pop), linewidth = 0.5) +
  #geom_smooth(aes(x.left, y, color = pop), method = lm, formula = y ~ poly(x, 8), se = F, linewidth = 0.5) +
  scale_x_continuous(breaks = xticks, labels = plot.marks(xticks), name = "Years ago", expand = c(0,0.02)) +
  scale_y_continuous(breaks = yticks, labels = plot.marks(yticks), name = expression(italic(N[e]))) +
  scale_color_manual(aesthetics = "color", values = colors, name = "Genogroup") +
  scale_fill_manual(values = colors, name = "95% CI", guide = "none") +
  theme_minimal() + 
  theme(text = element_text(color = "black", size = 7.5),
        aspect.ratio = 0.75,
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(colour = "black", size = 0.2),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        #legend.position = "none"
  )
out.plot.name <- paste0(path,"/population_plot_", mu_out, "_", gen, "_", random, ".tiff")
if (grepl("windows", os, ignore.case=T)){
  out.plot.name <- gsub("/", "\\\\", out.plot.name)
}
tiff(out.plot.name, units="cm", width = 8, height = 4, res = 600)
print(pop.plot)
dev.off()
out.plot.name.data <- paste0(path,"/population_plot_", mu_out, "_", gen, "_", random, ".Rdata")
save(final.table, file = out.plot.name.data)
