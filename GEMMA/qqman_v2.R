library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
os <- Sys.info()[['sysname']]
options(warn=-1)
s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Jul. 2022
Usage: Rscript qqman2.R -f ASSOC_FILE [-o OUTPUT_FILE_PREFIX] [-hl CHR:START-END] [-lm] [-h]\n
-f/--file: table file. Table seperates by tab with header line.
-lm/--linear: linear model.
-o/--output: output file name prefix without extention.
-hl/--highlight: highlight box for figure.
-h/--help: help.\n\n")
}
if (length(args) == 0){
  s.help()
  quit()
}
path <- c()
result_file <- c()
hl <- c()
lm <- c()
for (i in 1:length(args)){
  if (args[i] == '-f' || args[i] == '--file'){ #assoc file
    result_file <- args[i+1]
    if (!file.exists(result_file)){
      s.help()
      cat("-f: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-o' || args[i] == '--output'){ #output prefix name with path
    path <- args[i+1]
    name <- unlist(str_split(path, "/|\\\\"))
    check_path <- paste0(name[1:(length(name)-1)], collapse = "/", sep = "/")
    if (grepl("windows", os, ignore.case=T)){
      check_path <- gsub("/", "\\\\", check_path)
    }
    if (!dir.exists(check_path)){
      dir.create(check_path)
    }
  }
  if (args[i] == '-hl' || args[i] == '--highlight'){ #highlight box
    hl <- args[i+1]
  }
  if (args[i] == '-lm' || args[i] == '--linear'){ #highlight box
    lm <- 1
  }
  if (args[i] == '-h' || args[i] == '--help'){
    s.help()
    quit()
  }
  
}
if (length(result_file) == 0){
  s.help()
  quit()
}

if (length(path) == 0){
  path <- gsub(".assoc", "", result_file);
  if (grepl("windows", os, ignore.case=T)){
    path <- gsub("/", "\\\\", path)
  }
}

library(ggplot2)
library(gridExtra)

#calculate Bonferroni and BH threshold for GWAS (adopted from RAINBOWR package)
CalcThreshold <- function (input, sig.level = 0.05, method = "BH"){
  qvalue_tmp <- function(p) {
    smooth.df <- 3
    if (min(p) < 0 || max(p) > 1) {
      stop("P-values not in valid range.")
      return(0)
    }
    lambda <- seq(0, 0.9, 0.05)
    m <- length(p)
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    pi0 <- predict(spi0, x = max(lambda))$y
    pi0 <- min(pi0, 1)
    if (pi0 <= 0) {
      stop("The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
      idx <- sort.list(x)
      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl
      return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
      qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                          1)
    }
    return(qvalue)
  }
  input <- input[!is.na(input[, 4]), , drop = FALSE]
  input <- input[order(input[, 2], input[, 3]), ]
  method[!(method %in% c("BH", "Bonf"))] <- "BH"
  methods <- rep(method, each = length(sig.level))
  sig.levels <- rep(sig.level, length(method))
  n.thres <- length(methods)
  thresholds <- rep(NA, n.thres)
  for (thres.no in 1:n.thres) {
    method.now <- methods[thres.no]
    sig.level.now <- sig.levels[thres.no]
    if (method.now == "BH") {
      q.ans <- qvalue_tmp(10^(-input[, 4]))
      temp <- cbind(q.ans, input[, 4])
      temp <- temp[order(temp[, 1]), ]
      if (temp[1, 1] < sig.level.now) {
        temp2 <- tapply(temp[, 2], temp[, 1], mean)
        qvals <- as.numeric(rownames(temp2))
        x <- which.min(abs(qvals - sig.level.now))
        first <- max(1, x - 2)
        last <- min(x + 2, length(qvals))
        if ((last - first) < 4) {
          last <- first + 3
        }
        if (sum(is.na(qvals[first:last])) == 1) {
          qvals[last] <- mean(qvals[first + 1] + qvals[first + 2])
          temp2[last] <- mean(temp2[first + 1] + temp2[first + 2])
        }
        if (sum(is.na(qvals[first:last])) == 2) {
          qvals[(last - 1):last] <- quantile(qvals[first:(first + 1)], probs = c(1/3, 2/3))
          temp2[(last - 1):last] <- quantile(temp2[first:(first + 1)], probs = c(1/3, 2/3))
        }
        qvals <- sort(qvals)
        temp2 <- temp2[order(qvals)]
        splin <- smooth.spline(x = qvals[first:last], 
                               y = temp2[first:last], df = 3)
        threshold <- predict(splin, x = sig.level.now)$y
      } else {
        threshold <- NA
      }
    }
    if (method.now == "Bonf") {
      n.mark <- nrow(input)
      threshold <- -log10(sig.level.now/n.mark)
    }
    thresholds[thres.no] <- threshold
  }
  names(thresholds) <- paste0(methods, "_", sig.levels)
  return(thresholds)
}

#preparing data for manhattan plot
results <- read.table(result_file, header=T)
chr_len <- results %>% group_by(chr) %>% summarise(chr_len=max(pos, na.rm = TRUE))
chr_pos <- chr_len %>% mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>% select(-chr_len)
snp_pos <- chr_pos %>% left_join(results, ., by="chr") %>% arrange(chr, pos) %>% mutate(BPcum=pos+total)
caption.name <- unlist(str_split(path, "/|\\\\"))
X_axis <- snp_pos %>% group_by(chr) %>% summarize(center=(max(BPcum, na.rm = TRUE)+min(BPcum, na.rm = TRUE))/2)
max_y <- max(-log10(snp_pos$pval), na.rm = TRUE)
max_y <- max_y * 1.2


#preparing data for qq plot
ci <- 0.95
nSNPs <- nrow(snp_pos)
obs <- arrange(snp_pos, pval) 
obs <- -log10(obs$pval)
qq.data <- data.frame(
  observed = obs,
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
)
qq.data.sub <- qq.data %>% filter(expected <= 2) %>% sample_frac(0.01)
qq.data.sup <- qq.data %>% filter(expected > 2)
qq.data.small <- rbind(qq.data.sub, qq.data.sup)

#highlight boxes
box.begin <- c()
box.end <- c()
chr_pos <- as.data.frame(chr_pos)
if (length(hl) > 0){
  hl <- unlist(str_split(hl, ','))
  for (i in 1: length(hl)){
    chr <- unlist(str_split(hl[i], ':'))
    for (j in 1: nrow(chr_pos)){
      if (chr[1] == chr_pos[j,1]){
        box.current.begin <- as.numeric(chr_pos[j,2])
        box.current.end <- as.numeric(chr_pos[j,2])
      }
    }
    if (length(chr[2]) == 0){
      message("Format of highlighting box is wrong.")
      quit()
    }
    box.pos <- unlist(str_split(chr[2], '-'))
    if (length(box.pos[2]) == 0){
      message("Format of highlighting box is wrong.")
      quit()
    }
    box.current.begin <- box.current.begin + as.numeric(box.pos[1])
    box.current.end <- box.current.end + as.numeric(box.pos[2])
    box.begin <- c(box.begin, box.current.begin)
    box.end <- c(box.end, box.current.end)
  }
}

#two color theme
if (length(lm) == 0){
  threshold.data <- snp_pos[,c(2,1,3,14)]
} else {
  threshold.data <- snp_pos[,c(2,1,3,12)]
}
threshold.data$pval <- -log10(threshold.data$pval)
#BH <- as.numeric(CalcThreshold(threshold.data, sig.level = 0.05, method = "BH")) #red line
Bonf <- as.numeric(CalcThreshold(threshold.data, sig.level = 0.05, method = "Bonf")) #pink line
#Bonf <- -log10(0.05/nrow(snp_pos)) #pink line
man.plot <- ggplot()
if (length(hl) > 0){
  man.plot <- man.plot +  
    geom_rect(data = box.begin, mapping = aes(xmin = box.begin, xmax = box.end, ymin = max_y*0.98, ymax = max_y),
              alpha=1,
              #fill = "khaki1",
              fill = "grey30",
              inherit.aes = F)
}
man.plot <- man.plot +
  geom_point(data = snp_pos, mapping = aes(x=BPcum, y=-log10(pval), color=as.factor(chr)), size=1.5) + 
  scale_color_manual(values = rep(c("dodgerblue","gold"),length(chr_len$chr))) + 
  scale_x_continuous(label=X_axis$chr, breaks=X_axis$center) + 
  scale_y_continuous(limits=c(0,max_y), expand=c(0, 0)) +
  #geom_hline(yintercept=c(-log10(1e-5),BH,Bonf), color=c('darkred', 'red','pink1'),size = 0.5,linetype=c("longdash","longdash")) + 
  geom_hline(yintercept=c(-log10(1e-5),Bonf), color=c('red','pink1'),size = 0.5,linetype=c("longdash","longdash")) + 
  guides(color = "none") + theme_classic() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.caption = element_text(face = "italic"),
        text = element_text(color = "black", face = "bold", size = 18)) +
  labs(x = "Chromosome", y = expression(bold(paste("-log"["10"],"(", bolditalic(p),")"))), caption = caption.name[length(caption.name)])

qq.plot <- ggplot(qq.data.small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "blue", alpha = 0.25) +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "blue", lineend = "round") +
  geom_point(color = "black", size = 1.5, alpha = 0.8) +
  labs(x = expression(bold(paste("Expected -log"["10"],"(", bolditalic(p),")"))),
       y = expression(bold(paste("Observed -log"["10"],"(", bolditalic(p),")")))) +
  theme_minimal() + theme(aspect.ratio = 1, 
                          panel.border = element_rect(color = "black", fill = NA, size = 1.7), 
                          panel.grid = element_blank(),
                          axis.ticks = element_line(colour = "black", size = 1.5),
                          axis.title = element_text(face = "bold"),
                          axis.text = element_text(color = "black"),
                          legend.title = element_text(face = "bold"),
                          text = element_text(color = "black", face = "bold", size = 18))
tiff(paste0(path,".qqman_plot.tiff"), units = "in", pointsize = 12, res = 300, bg = "white", compression = c("none"), width = 16, height = 4)
grid.arrange(qq.plot, man.plot, nrow = 1, widths = c(3.5,12), heights = 4)
dev.off()

