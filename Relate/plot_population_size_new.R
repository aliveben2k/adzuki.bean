
#check if packages are installed
#list.of.packages <- c("ggplot2", "grid", "gridExtra", "dplyr")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(gridExtra)
library(dplyr)
os <- Sys.info()[['sysname']]

#args <- commandArgs(trailingOnly = T)
#colors <- c("#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
colors <- c("#4b8bcb","gold","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")

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

#define seperation time for RCCR
det.range <- function(col.x = NULL, col.y = NULL, point = 0.5, side = "left") {
  minus <- length(col.x)
  max.y <- max(col.y, na.rm = T)
  meet.one <- 0
  sep.side.x <- c()
  sep.side.y <- c()
  for (a in 1:length(col.x)){
    if (col.y[minus] == max.y){
      meet.one <- 1
    }
    if (col.y[minus] < point && meet.one == 1){
      if (side == "right"){
        sep.side.x <- col.x[minus+1]
        sep.side.y <- col.y[minus+1]
      } else {
        sep.side.x <- col.x[minus]
        sep.side.y <- col.y[minus]        
      }
      break
    }
    minus <- minus-1
  }
  return(c(sep.side.x, sep.side.y))
}

args <- commandArgs(trailingOnly = TRUE)
dirs <- list.dirs(args[1], full.names = T, recursive = T)
years_per_gen <- as.numeric(args[2])
random <- args[3]
files <- c()
for (k in 1:length(dirs)){
  files.tmp <- list.files(dirs[k], pattern="pairwise\\.coal$", full.names=TRUE)
  files.tmp <- files.tmp[!grepl("combined", files.tmp)]
  tmp <- unlist(strsplit(files.tmp[1],split="/"))
  folder.name <- tmp[length(tmp)-1]
  files <- c(files, files.tmp[!grepl(paste0(folder.name, "_"), files.tmp)])
}
if (length(files) == 0){
  quit("no")
}

#combine all data
gens.all <- c()
pairs <- list()
groups <- c()
for (i in 1:length(files)){
  groups <- read.table(file = files[i], header = F, nrows = 1) #get population name and order
  gens <- read.table(file = files[i], header = F, nrows = 1, skip = 1) #get generation defined in the files
  values <- read.table(file = files[i], header = F, skip = 2)
  gens.all <- rbind(gens.all, gens)
  list.idx <- 1
  for (pop1_idx in 0:(ncol(groups)-1)){
    for (pop2_idx in pop1_idx:(ncol(groups)-1)){
      tmp.value <- values %>% filter(V1 == pop1_idx & V2 == pop2_idx)
      if (i == 1){
        pairs[[list.idx]] <- as.data.frame(tmp.value)
      } else {
        pairs[[list.idx]] <- rbind(pairs[[list.idx]], as.data.frame(tmp.value))
      }
      list.idx <- list.idx + 1
    }
  }
  time <- t(years_per_gen*gens)
  colnames(time) <- "time"
}
#averaging
combined.table <- c()
ci.upper.table <- c()
ci.lower.table <- c()
for (i in 1:length(pairs)){
  avg.values <- colMeans(pairs[[i]])
  combined.table <- rbind(combined.table, avg.values)
  #calculate 95% CI
  sd.values <- sqrt(diag(var(pairs[[i]])))
  se.values <- sd.values/sqrt(nrow(pairs[[i]]))
  alpha <- 0.05
  if (length(files) > 1){
    degrees.freedom <- nrow(pairs[[i]]) - 1
    t.score <- qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
    margin.error <- t.score * se.values
  } else {
    margin.error <- 0
  }
  upper.bound <- avg.values + margin.error
  lower.bound <- avg.values - margin.error
  upper.bound[c(1,2)] <- avg.values[c(1,2)]
  ci.upper.table <- rbind(ci.upper.table, upper.bound)
  lower.bound[c(1,2)] <- avg.values[c(1,2)]
  ci.lower.table <- rbind(ci.lower.table, lower.bound)
}
gens.all <- data.frame(colMeans(gens.all))

#read in population size
time <- years_per_gen*gens.all
colnames(time) <- "time"
pop_size <- data.frame(time = numeric(0), pop_size = numeric(0), groups = numeric(0))
pop_size.ci.upper <- data.frame(time = numeric(0), upper = numeric(0), groups = numeric(0))
pop_size.ci.lower <- data.frame(time = numeric(0), lower = numeric(0), groups = numeric(0))
rccr <- data.frame(time = numeric(0), coal = numeric(0), groups = numeric(0))
rccr.ci.upper <- data.frame(time = numeric(0), upper = numeric(0), groups = numeric(0))
rccr.ci.lower <- data.frame(time = numeric(0), lower = numeric(0), groups = numeric(0))
#ccr <- data.frame(time = numeric(0), ccr = numeric(0), groups = numeric(0))
num_pops <- length(groups)

for(p1 in 1:num_pops){
  for(p2 in p1:num_pops){
    c <- as.data.frame(combined.table) %>% filter(V1 == (p1-1) & V2 == (p2-1)) %>% select(-c(1,2))
    c.upper <- as.data.frame(ci.upper.table) %>% filter(V1 == (p1-1) & V2 == (p2-1)) %>% select(-c(1,2))
    c.lower <- as.data.frame(ci.lower.table) %>% filter(V1 == (p1-1) & V2 == (p2-1)) %>% select(-c(1,2))
    c <- t(c)
    c.upper <- t(c.upper)
    c.lower <- t(c.lower)
    colnames(c) <- "pop_size"
    colnames(c.upper) <- "upper"
    colnames(c.lower) <- "lower"
    str <- rep(paste(groups[p1]),length(c))
    if (groups[p1] == groups[p2]){
      pop_size <- rbind(pop_size, data.frame(time, 0.5/c, str))
      pop_size.ci.upper <- rbind(pop_size.ci.upper, data.frame(time, 0.5/c.upper, str))
      pop_size.ci.lower <- rbind(pop_size.ci.lower, data.frame(time, 0.5/c.lower, str))
    }
    rccr <- rbind(rccr, data.frame(time, c, paste(groups[p1]," - ",groups[p2], sep = "")))
    rccr.ci.upper <- rbind(rccr.ci.upper, data.frame(time, c.upper, paste(groups[p1]," - ",groups[p2], sep = "")))
    rccr.ci.lower <- rbind(rccr.ci.lower, data.frame(time, c.lower, paste(groups[p1]," - ",groups[p2], sep = "")))
    #str_ccr <- rep(paste(groups[p1]," - ",groups[p2], sep = ""),length(c))
    #ccr <- rbind(ccr, data.frame(time = time, ccr = c, groups = str_ccr))
  }
}
colnames(pop_size)[3] <- "groups"
colnames(pop_size.ci.upper)[3] <- "groups"
colnames(pop_size.ci.lower)[3] <- "groups"
colnames(rccr)[3] <- "groups"
colnames(rccr.ci.upper)[3] <- "groups"
colnames(rccr.ci.lower)[3] <- "groups"
rccr <- data.frame(cbind(time = rccr$time, coal = rccr$pop_size, upper = rccr.ci.upper$upper, lower = rccr.ci.lower$lower, groups = rccr$groups))

#define number of colors in the plot
dis.cats <- length(groups)
if (dis.cats >= 11){
  if(!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  colors <- getPalette(dis.cats)
} else {
  colors <- colors[1:dis.cats]
}

#plot
pop_size <- pop_size[which(pop_size$pop_size != "NaN"),]
pop_size.ci.upper <- pop_size.ci.upper[which(pop_size.ci.upper$upper != "NaN"),]
pop_size.ci.lower <- pop_size.ci.lower[which(pop_size.ci.lower$lower != "NaN"),]
#Va use
# pop_size <- pop_size[!grepl("SOUTH_LR|CN_LRS",pop_size$groups),]
# pop_size.ci.upper <- pop_size.ci.upper[!grepl("SOUTH_LR|CN_LRS",pop_size.ci.upper$groups),]
# pop_size.ci.lower <- pop_size.ci.lower[!grepl("SOUTH_LR|CN_LRS",pop_size.ci.lower$groups),]
# pop_size$groups <- factor(pop_size$groups, levels = c("JP_LR", "CN_LRN","JP_WL", "CN_WL", "SOUTH_WL2", "SOUTH_WL1"))
#end of Va use
tmax <- log10(pop_size$time[which(pop_size$time < Inf)])
x.min.ticks <- sort(tmax[which(tmax > -Inf)])[1]
xticks <- plot.ticks(x.min.ticks, max(tmax, na.rm = T))
yticks <- plot.ticks(0, log10(max(pop_size.ci.upper$upper, na.rm = T)))
pop_size$time.end <- c(unlist(pop_size$time[2:nrow(pop_size)]), "Inf")
pop_size$time.end[which(pop_size$time.end == 0)] <- "Inf"
pop_size <- cbind(pop_size, pop_size.ci.upper$upper, pop_size.ci.lower$lower)
colnames(pop_size)[c((ncol(pop_size)-1),ncol(pop_size))] <- c("lower","upper")
out.ps.name <- paste0(dirs,"/", random, "_pop_size.txt")
write.table(pop_size, file = out.ps.name, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)
options(scipen=999)
pop.plot <- ggplot(pop_size) +
  geom_rect(aes(xmin = log10(as.numeric(time)), xmax = log10(as.numeric(time.end)), ymin = log10(as.numeric(lower)), ymax = log10(as.numeric(upper)), fill = groups), alpha = 0.5) +
  geom_step(aes(log10(time), log10(pop_size), color = groups), linewidth = 0.3) +
  #geom_point(size = 0.5, alpha = 0.25) +
  #geom_smooth(method = lm, formula = y ~ poly(x, 8), se = F, linewidth = 0.5) +
  scale_x_continuous(breaks = xticks, labels = plot.marks(xticks), name = "Years ago", limits = c(3,5.8)) +
  #scale_x_continuous(breaks = xticks, labels = plot.marks(xticks), name = "Years ago") +
  scale_y_continuous(breaks = yticks, labels = plot.marks(yticks), name = expression(italic(N[e]))) +
  scale_color_manual(values = colors, name = "Genogroup") +
  scale_fill_manual(values = colors, name = "95% CI", guide = "none") +
  theme_minimal() +
  theme(aspect.ratio = 0.75,
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.line = element_line(colour = "black", size = 0.2),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(color = "black", size = 7.5),
        legend.position = "right",
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(color = "black", size = 7.5))

out.ps.plot.name <- paste0(dirs,"/", random, "_popsize_plot.tiff")
tiff(out.ps.plot.name, units="cm", width = 8, height = 4, res = 600)
print(pop.plot)
dev.off()

rccr.plots <- list()
range.box <- list()
rccr.tables <- list()
k <- 1 #count plot numbers
for (i in 1:num_pops){
  ##Va use
  # if (groups[i] == "SOUTH_LR" || groups[i] == "CN_LRS"){
  #   next
  # }
  # for (j in 2:(num_pops-1)){
  #   if (groups[j] == "SOUTH_LR" || groups[j] == "CN_LRS"){
  #     next
  #   }
  ##end of Va use
    if (groups[i] == groups[j]){
      next
    }
    rccr.table <- c()
    pop1_ccr <- rccr[which(rccr$groups == paste0(groups[i], " - ", groups[i])),]
    pop1_ccr <- pop1_ccr %>% filter(coal != "NaN") %>% mutate_at(1:(ncol(pop1_ccr)-1), as.numeric)
    pop2_ccr <- rccr[which(rccr$groups == paste0(groups[j], " - ", groups[j])),]
    pop2_ccr <- pop2_ccr %>% filter(coal != "NaN") %>% mutate_at(1:(ncol(pop2_ccr)-1), as.numeric)
    pop12_ccr <- rccr[which(rccr$groups == paste0(groups[i], " - ", groups[j])),]
    pop12_ccr <- pop12_ccr %>% filter(coal != "NaN") %>% mutate_at(1:(ncol(pop12_ccr)-1), as.numeric)
    if (nrow(pop12_ccr) == 0){
      pop12_ccr <- rccr[which(rccr$groups == paste0(groups[j], " - ", groups[i])),]
      pop12_ccr <- pop12_ccr %>% filter(coal != "NaN") %>% mutate_at(1:(ncol(pop12_ccr)-1), as.numeric)
    }
    time.row <- min(c(nrow(pop1_ccr),nrow(pop12_ccr),nrow(pop2_ccr)))
    pop1_ccr <- pop1_ccr[1:time.row,]
    pop2_ccr <- pop2_ccr[1:time.row,]
    pop12_ccr <- pop12_ccr[1:time.row,]
    rccr.table <- data.frame(cbind(pop1_ccr[,1:2],pop12_ccr[,2],pop2_ccr[,2],pop1_ccr[,3],pop12_ccr[,3],pop2_ccr[,3],pop1_ccr[,4],pop12_ccr[,4],pop2_ccr[,4]))
    colnames(rccr.table)[2:ncol(rccr.table)] <- c("lambda_00","lambda_01","lambda_11","lambda_00_low","lambda_01_low","lambda_11_low","lambda_00_up","lambda_01_up","lambda_11_up")
    x <- log10(as.numeric(rccr.table$time))
    y <- (2*rccr.table$lambda_01)/(rccr.table$lambda_00+rccr.table$lambda_11)
    y.lower <- (2*rccr.table$lambda_01_low)/(rccr.table$lambda_00_low+rccr.table$lambda_11_low)
    y.upper <- (2*rccr.table$lambda_01_up)/(rccr.table$lambda_00_up+rccr.table$lambda_11_up)
    rccr.table <- cbind(rccr.table, x, y, y.lower, y.upper)
    tmax <- rccr.table$x[which(rccr.table$x < Inf)]
    x.min.ticks <- sort(tmax[which(tmax > -Inf)])[1]
    xticks <- plot.ticks(x.min.ticks,max(tmax, na.rm = T))
    options(scipen=999)
    x.max <- c()
    #optional
    ##define the right boundary of the x-axis for each population
    # all.names <- c("CN_LRN","CN_WL","JP_LR","JP_WL","SOUTH_WL1","SOUTH_WL2")
    # all.limit <- c(30000, 200000, 30000, 70000, 1000000, 1000000)
    # all.limit <- log10(all.limit)
    # for (z in 1:length(all.names)){
    #   if (groups[i] == all.names[z]){
    #     x.max.1 <- all.limit[z]
    #   }
    #   if (groups[j] == all.names[z]){
    #     x.max.2 <- all.limit[z]
    #   }
    # }
    # x.max <- min(x.max.1, x.max.2)
    # rccr.table <- rccr.table[rccr.table$x <= x.max,]
    ##end of optional code
    if (length(x.max) == 0){
      x.max <- rccr.table$x[nrow(rccr.table)]
    }
    scale.y <- max(rccr.table$y, na.rm = T)
    rccr.table$y <- rccr.table$y/scale.y #scaling
    rccr.table$y.lower <- rccr.table$y.lower/scale.y
    rccr.table$y.upper <- rccr.table$y.upper/scale.y
    out.rccr.name <- paste0(dirs,"/", random, '_', groups[i], '-', groups[j],"_rccr.txt")
    write.table(rccr.table, file = out.rccr.name, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)
    rccr.tables[[k]] <- rccr.table
    sep.point <- c()
    sep.start <- c()
    mid.point <- 0.5
    mid.start <- 0.8
    sep.start <- det.range(col.x = rccr.table$x, col.y = rccr.table$y, point = mid.start, side = "right")
    sep.point <- det.range(col.x = rccr.table$x, col.y = rccr.table$y, point = mid.point, side = "left")
    sep.start.ci.upper <- det.range(col.x = rccr.table$x, col.y = rccr.table$y.upper, point = mid.start, side = "right")
    sep.start.ci.lower <- det.range(col.x = rccr.table$x, col.y = rccr.table$y.lower, point = mid.start, side = "right")
    sep.point.ci.upper <- det.range(col.x = rccr.table$x, col.y = rccr.table$y.upper, point = mid.point, side = "left")
    sep.point.ci.lower <- det.range(col.x = rccr.table$x, col.y = rccr.table$y.lower, point = mid.point, side = "left")
    sep.points <- approx(c(sep.point[2], sep.start[2]), c(sep.point[1], sep.start[1]), xout = c(mid.point,mid.start))
    sep.points.upper <- approx(c(sep.point.ci.upper[2], sep.start.ci.upper[2]), c(sep.point.ci.upper[1], sep.start.ci.upper[1]), xout = c(mid.point,mid.start))
    sep.points.lower <- approx(c(sep.point.ci.lower[2], sep.start.ci.lower[2]), c(sep.point.ci.lower[1], sep.start.ci.lower[1]), xout = c(mid.point,mid.start))
    range.box[[k]] <- as.data.frame(c(sep.points[["y"]][1], sep.points[["y"]][2]))
    range.box.ci <- as.data.frame(cbind(c(sep.points.upper[["y"]][1], sep.points.upper[["y"]][2]),c(sep.points.lower[["y"]][1], sep.points.lower[["y"]][2])))
    pop_pair <- paste0(groups[i], '_', groups[j])
    range.box[[k]] <- cbind(range.box[[k]], range.box.ci, rep(pop_pair, nrow(range.box[[k]])), c(mid.point,mid.start))
    colnames(range.box[[k]]) <- c("time.point","upper","lower","pair", "rccr")
    rownames(range.box[[k]]) <- c(mid.point,mid.start)
    sep.plot <- ggplot() +
      geom_rect(data = range.box[[k]], xmin = range.box[[k]][1,2], xmax = range.box[[k]][1,3], ymin = 0, ymax = Inf, fill = "#f6ba75", alpha = 0.7) + # rccr = 0.5, 95% CI
      geom_rect(data = range.box[[k]], xmin = range.box[[k]][2,2], xmax = range.box[[k]][2,3], ymin = 0, ymax = Inf, fill = "#f6ba75", alpha = 0.7) + #rccr = 0.8, 95% CI
      geom_rect(data = range.box[[k]], xmin = range.box[[k]][1,1], xmax = range.box[[k]][2,1], ymin = 0, ymax = Inf, fill = "#fff2a7")
    sep.plot <- sep.plot +
      #geom_rect(data = rccr.table, aes(xmin = mean.x.left, xmax = mean.x.right, ymin = lower.bound, ymax = upper.bound), fill = "#4b8bcb", alpha = 0.7)
      geom_ribbon(data = rccr.table, aes(x = x, y = y, xmin = x, xmax = x, ymin = y.lower, ymax = y.upper), fill = "#4b8bcb", color = NA, alpha = 0.7, show.legend = NA)
    sep.plot <- sep.plot + 
      geom_line(data = rccr.table, aes(x, y), color = "black", size = 0.5)
    #geom_smooth(rccr.table, mapping = aes(x, y), color = "black", method = lm, formula = y ~ poly(x, 8), se = F, linewidth = 1)
    sep.plot <- sep.plot + 
      scale_x_continuous(breaks = xticks, 
                         labels = plot.marks(xticks),
                         #limits = c(2.95,x.max), #optional code
                         expand = c(0,0),
                         name = "Years ago") +
      scale_y_continuous(name = "Relative CCR",
                         expand = c(0,0),
                         breaks=seq(0,2,0.25)) +
      labs(title=paste(groups[i], 'vs.', groups[j], collapse = " ")) +
      theme_minimal() +
      theme(aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid = element_blank(),
            axis.ticks = element_line(colour = "black", size = 0.2),
            axis.text = element_text(color = "black", size = 7.5),
            legend.position = "none",
            text = element_text(color = "black", size = 7.5))
    options(scipen=1)
    rccr.plots[[k]] <- sep.plot
    k <- k+1
  # } #Va use
}

#for outputting separation times
range.box.all <- c()
for(i in 1:length(range.box)){
  range.box.all <- rbind(range.box.all, range.box[[i]])
}
range.box.all$time.point <- 10^range.box.all$time.point
range.box.all$upper <- 10^range.box.all$upper
range.box.all$lower <- 10^range.box.all$lower
out.box.name <- paste0(dirs,"/", random,"_popsize.pairwise.sep_range.txt")
write.table(range.box.all, file = out.box.name, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)


plot.heights = sum(rep(4.5, ceiling(k/3)))
#out.plot.name <- paste(filename,".rccr.tiff",sep = "")
out.plot.name <- paste0(dirs,"/", random,"_popsize.pairwise.rccr.tiff")
if (grepl("windows", os, ignore.case=T)){
  out.plot.name <- gsub("/", "\\\\", out.plot.name)
}
tiff(out.plot.name, width=18, height=plot.heights, res = 600, units = "cm")
grid.arrange(grobs = rccr.plots,
             ncol = 3,
             widths = rep(6, 3), 
             heights = rep(4.5, ceiling(k/3)))
dev.off()
save(pop_size, pop_size.ci.lower, pop_size.ci.upper, rccr.tables, range.box.all, file = paste0(dirs,"/", random,"_popsize.analysis.all.Rdata"))


