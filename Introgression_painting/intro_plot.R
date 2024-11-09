library(ggplot2)
library(dplyr)

s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Oct. 2024
Usage: Rscript intro_plot.R -p PATH -t FILE -gi FILE [-ci] [-p1c COLOR] [-p2c COLOR]\n
-p/--path: path of rda files generated from new_intro_count.R.
-t/--trios: trios information file. (samples are seperate by tab)
  Format:
  Parent_Pop1 Sample1 Sample2...
  Parent_Pop2 Sample3 Sample4...
  Test_Pop Sample5 Sample6...
-gi/--genome_info: genome information generated from vcf2trios_thread.pl (*genome_info.txt).
-ci/--ci: show 95% confidence interval. Default: False.
-p1c/--p1_color: the color to indicate the ratio from ancestor 1. Default: blue.
-p2c/--p2_color: the color to indicate the ratio from ancestor 2. Default: yellow.\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
path <- c()
pop.info <- c()
geno.info.file <- c()
CI.switch <- 0
p1.color = "#4b8bcb"
p2.color = "gold"
for (i in 1:length(args)){
  if (args[i] == '-t' || args[i] == '--trios'){ #assume the first two rows are two parents, and the third row is the test population
    pop.info <- as.character(args[i+1])
    if (!file.exists(pop.info)){
      cat("-t: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-gi' || args[i] == '--genome_info'){ #genome information file containing chromosome names and lengths
    geno.info.file <- as.character(args[i+1])
    if (!file.exists(geno.info.file)){
      cat("-gi: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-p' || args[i] == '--path'){ #path of rda files
    path <- as.character(args[i+1])
    if (!file.exists(path)){
      cat("-p: path does not exist.\n")
      quit()
    }
    if (grepl("/$", path)){
      path <- sub("/$", "", path)
    }
  }
  if (args[i] == '-ci' || args[i] == '--ci'){
    CI.switch <- 1
  }
  if (args[i] == '-p1c' || args[i] == '--p1_color'){
    p1.color <- as.character(args[i+1])
  }
  if (args[i] == '-p2c' || args[i] == '--p2_color'){
    p2.color <- as.character(args[i+1])
  }
}

dirs <- list.dirs(path, full.names = T, recursive = T)
files <- c()
for (k in 1:length(dirs)){
  files <- list.files(dirs[k], pattern="introgression\\.rda$", full.names=TRUE)
}
if (length(files) == 0){
  quit("no")
}

#read the genome infomation
geno.info <- read.table(geno.info.file, header = T, sep = "\t")
geno.info <- as.data.frame(geno.info)

#read the population info
lines <- readLines(pop.info)
trio <- c()
for (i in 1:length(lines)){
  line.elements <- unlist(strsplit(lines[i],"\t"))
  trio <- c(trio, line.elements[1]) #the file will contain trio[3] (child) information
}

#load all file contents into a list
all.samples <- list()
curr.sample <- c()
for (i in 1:length(files)){
  load(files[i])
  curr.sample <- colnames(curr.unique.win)[3]
  if (!is.null(all.samples[[curr.sample]])){
    all.samples[[curr.sample]] <- rbind(all.samples[[curr.sample]], curr.unique.win)
  } else {
    all.samples[[curr.sample]] <- curr.unique.win
  }
}

for (sample in names(all.samples)){
  curr.unique.win <- all.samples[[sample]]
  #sort by chromosome order
  curr.unique.win <- curr.unique.win %>% mutate(order = match(Chr, geno.info$Chr)) %>%
    arrange(order) %>%
    select(-order)
  all.chrs <- unique(curr.unique.win$Chr)
  geno.info <- geno.info[geno.info$Chr %in% all.chrs,]
  chr.pos <- geno.info %>% mutate(total = cumsum(as.numeric(Length)) - Length) %>% select(-Length)
  curr.unique.win <- chr.pos %>% left_join(curr.unique.win, ., by="Chr") %>% arrange(Chr, Pos) %>% mutate(BPcum=Pos+total)
  X_axis <- curr.unique.win %>% group_by(Chr) %>% summarize(center=(max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE))/2)
  X_axis$Chr <- gsub("chr|chr0", "", X_axis$Chr, ignore.case = T)
  X_lines <- chr.pos[2:nrow(chr.pos),2]
  curr.plot <- ggplot(curr.unique.win, aes(x = BPcum, y = curr.unique.win[,3])) +
    geom_ribbon(aes(ymin = 0, ymax = curr.unique.win[,3]), fill = p1.color, alpha = 1) +
    geom_ribbon(aes(ymin = curr.unique.win[,3], ymax = 1), fill = p2.color, alpha = 1) +
    geom_vline(xintercept = X_lines, color = "white") +
    scale_x_continuous(label=X_axis$Chr, breaks=X_axis$center, expand = c(0.02, 0)) +
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks = c(0,0.5,1)) +
    labs(x = "Chromosome", y = "Ratio", title = colnames(curr.unique.win)[3]) +
    guides(color = "none") + theme_minimal() + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          axis.text = element_text(color = "black", size = 7.5),
          axis.line.y = element_line(colour = "black", size = 0.2),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.2),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.ticks.x = element_blank(),
          #axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0,vjust = 0.5),
          #axis.text.x = element_blank(),
          plot.caption = element_text(face = "italic"),
          plot.title = element_text(color = "black", size = 7.5),
          text = element_text(color = "black", size = 7.5))
  if (CI.switch == 1){
    curr.plot <- curr.plot +
      geom_ribbon(aes(ymin = curr.unique.win[,4], ymax = curr.unique.win[,5]), fill = "white", alpha = 0.25)
  }  
  tiff(paste0(path,"/", colnames(curr.unique.win)[3],"_introgression.tiff"), units = "cm",res = 600, width = 15, height = 2)
  print(curr.plot)
  dev.off()
}
