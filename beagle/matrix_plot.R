#args[1]: matrix file
#args[2]: pop file
library(reshape2)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
dat <- read.table(args[1], sep="\t", header=TRUE)
melted_dat <- melt(dat, id.vars = "X", variable.name = "Y")
list <- read.table(args[2], sep="\t", header=F)
min <- min(melted_dat$value)
max <- max(melted_dat$value)
avg <- (min+max)/2
names <- dat[1:length(dat),1]
#melted_dat$X <- factor(melted_dat$X, levels = names)
melted_dat$X <- factor(melted_dat$X, levels = list$V1)
melted_dat$Y <- factor(melted_dat$Y, levels = rev(list$V1))
pop_tmp <- list$V2[1]
position.break <- c()
vline.break <- c()
pop.name <- c()
for (i in 1:nrow(list)){
  if (i == 1){
    start.pos <- 1
    pop.name <- pop_tmp
  }
  if (list$V2[i] != pop_tmp){
    midpoint <- as.integer((start.pos+i-1)/2 + start.pos - 1)
    position.break <- c(position.break, midpoint)
    pop_tmp <- list$V2[i]
    pop.name <- c(pop.name, pop_tmp)
    start.pos <- i
    vline.break <- c(vline.break, i-0.5)
  }
  if (i == nrow(list)){
    midpoint <- as.integer((start.pos+nrow(list))/2 + start.pos - 1)
    position.break <- c(position.break, midpoint)
  }
}
hline.break <- nrow(list) - vline.break
hline.break <- sort(hline.break)
hline.break <- hline.break + 1
tiff(paste0(args[1],".tiff"), units = "in", res = 300, bg = "white", compression = c("none"), width = 17, height = 15)
ggplot(melted_dat, aes(x=X, y=Y, fill=value)) + 
  geom_tile() + 
  geom_vline(xintercept=vline.break, size = 1, color = "white", linetype="dotted") +
  geom_hline(yintercept=hline.break, size = 1, color = "white", linetype="dotted") +
  scale_fill_gradient2(name = "Total cM length", 
                      low = "navy", 
                      high = "moccasin",
                      mid = "gold",
                      midpoint = 10,
                      limits = c(0, 20), 
                      na.value = "moccasin") +
  labs(x="", y="") +
  guides(fill = guide_colourbar(title.vjust = 4,label.theme = element_text(colour = "black", size = 36))) +
  scale_x_discrete(labels = pop.name, breaks = position.break, position = "top") +
  #scale_x_continuous(labels = pop.name, breaks = position.break, expand = c(0,0), position = "top") +
  scale_y_discrete(labels = rev(pop.name), breaks = position.break, position = "right") +
  theme(aspect.ratio = 1,
        legend.key.size = unit(1, 'in'),
        legend.title = element_text(size = 38, face = "bold", color = "black"),
        text = element_text(color = "black", size = 36),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, color = "black", size = 34, face = "bold"),
        axis.text.y = element_text(size = 34, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())
dev.off()

rownames(dat) <- dat[,1]
dat <- dat[,2:ncol(dat)]
#list <- read.table(args[2], sep="\t", header=F)
pop <- unique(list$V2)
pop.values <- c()
pop.matrix.name <- c()
for (k in 1:nrow(list)){ #list
  for (l in 1:nrow(dat)){
    if (rownames(dat)[l] == list[k,1]){
      pop.idx <- grep(list[k,2], pop)
      pop.matrix.name <- c(pop.matrix.name, pop.idx)
    }
  }  
}
dat <- cbind(dat, pop.matrix.name)
pop.matrix.name <- c(pop.matrix.name, NA)
dat <- rbind(dat, pop.matrix.name)

pop.matrix <- rep(NA, length(pop)*length(pop))
pop.matrix <- matrix(pop.matrix, nrow=length(pop), ncol=length(pop), byrow = T)
pop.matrix.cnt <- rep(NA, length(pop)*length(pop))
pop.matrix.cnt <- matrix(pop.matrix.cnt, nrow=length(pop), ncol=length(pop), byrow = T)
for (i in 1:(nrow(dat)-1)){ #matrix row
  for (j in 1:(ncol(dat)-1)){ #matrix col
      pop.matrix.row <- dat[i, ncol(dat)]
      pop.matrix.col <- dat[nrow(dat), j]
      if (is.na(pop.matrix[pop.matrix.row,pop.matrix.col])){
        pop.matrix[pop.matrix.row,pop.matrix.col] <- 0
        pop.matrix.cnt[pop.matrix.row,pop.matrix.col] <- 0
      }
      pop.matrix[pop.matrix.row,pop.matrix.col] <- pop.matrix[pop.matrix.row,pop.matrix.col] + dat[i,j]
      pop.matrix.cnt[pop.matrix.row,pop.matrix.col] <- pop.matrix.cnt[pop.matrix.row,pop.matrix.col] + 1
  }
}

for (i in 1:nrow(pop.matrix)){
  for (j in 1:ncol(pop.matrix)){
    if (!is.na(pop.matrix[i,j])){
      pop.matrix[i,j] <- pop.matrix[i,j]/pop.matrix.cnt[i,j]
    }
  }
}
pop.matrix <- cbind(pop, pop.matrix)
colnames(pop.matrix) <- c("X", pop)
rownames(pop.matrix) <- pop
pop.matrix <- pop.matrix[,2:ncol(pop.matrix)]
melted_dat <- melt(pop.matrix, id.vars = 0, variable.name = "Y")
#melted_dat <- melted_dat[(length(pop)+1):nrow(melted_dat),]
melted_dat$value <- as.numeric(melted_dat$value)
max <- max(melted_dat$value, na.rm = T)
min <- min(melted_dat$value, na.rm = T)
avg <- (min+max)/2
melted_dat$Var1 <- factor(melted_dat$Var1, levels = pop)

tiff(paste0(args[1],".pop.tiff"), units = "in", res = 300, bg = "white", compression = c("none"), width = 8.5, height = 7.5)
ggplot(melted_dat, aes(x=Var1, y=rev(Var2), fill=value)) + 
  geom_tile() +
  #scale_x_continuous(labels = pop.name, breaks = position.break, expand = c(0,0), position = "top") +
  scale_fill_gradient2(name = "Total cM length", 
                       low = "navy", 
                       high = "moccasin",
                       mid = "gold",
                       midpoint = 10,
                       limits = c(0, 20), 
                       na.value = "moccasin") +
  scale_x_discrete(labels = pop.name) +
  scale_y_discrete(labels = rev(pop.name)) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.5, 'in'),
        legend.title = element_text(size = 20, face = "bold", color = "black"),
        text = element_text(color = "black", size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 16, face = "bold", margin = margin(b = 5)),
        axis.text.y = element_text(size = 16, hjust = 1, face = "bold", color = "black", margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())
dev.off()