library(ggplot2)
library(dplyr)

s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Oct. 2024
Usage: Rscript new_intro_count.R -g FILE -t FILE -gi FILE [-p PATH] [-m NUM] [-d NUM] [-w NUM] [-s NUM]\n
-g/--genotype: genotype file from vcf2trios_thread.pl (*.trios.gz).
-t/--trios: trios information file. (samples are seperate by tab)
  Format:
  Parent_Pop1 Sample1 Sample2...
  Parent_Pop2 Sample3 Sample4...
  Test_Pop Sample5 Sample6...
-gi/--genome_info: genome information generated from vcf2trios_thread.pl (*genome_info.txt).
-p/--path: output path.
-m/--missing: missing rate threshold (0~1). Default: 0.8
-d/--diff_threshold: ratio of genotype difference between two parents (0-1). Default: 0.8
-w/--window: window size for plot (unit: SNP number). Default: 1000
-s/--step: window size for plot (unit: SNP number). Default: 500\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
geno.file <- c()
pop.info <- c()
geno.info.file <- c()
path <- "."
window.size <- 1000 #SNP number
step.size <- 500 #overlapping SNP number
missing.rate <- 0.2 #missing rate threshold
threshold.diff <- 0.8 #genotype difference between p1 and p2
for (i in 1:length(args)){
  if (args[i] == '-g' || args[i] == '--genotype'){ #genotype file from vcf2trios_thread.pl
    geno.file <- as.character(args[i+1])
    if (!file.exists(geno.file)){
      cat("-g: file does not exist.\n")
      quit()
    }
  }
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
  if (args[i] == '-p' || args[i] == '--path'){ #path for output files
    path <- as.character(args[i+1])
    if (!file.exists(path)){
      cat("-p: path does not exist.\n")
      quit()
    }
    if (grepl("/$", path)){
      path <- sub("/$", "", path)
    }
  }
  if (args[i] == '-m' || args[i] == '--missing'){
    missing.rate <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", missing.rate) || missing.rate > 1 || missing.rate < 0){
      cat("-m: only a number between 0~1 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-d' || args[i] == '--diff_threshold'){
    threshold.diff <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", threshold.diff) || threshold.diff > 1 || threshold.diff < 0){
      cat("-d: the threshold between 0~1 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-w' || args[i] == '--window'){
    window.size <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", window.size) || window.size < 10){
      cat("-w: only numbers >= 10 is allowed.\n")
      quit()
    }
  }
  if (args[i] == '-s' || args[i] == '--step'){
    step.size <- as.numeric(args[i+1])
    if (grepl("[^0-9.]", step.size) || step.size < 10){
      cat("-s: only numbers >= 10 is allowed.\n")
      quit()
    }
  }
}

# check window size and step size
if (step.size >= window.size){
    cat("-s: the step size cannot be equal or larger than the window size.\n")
    quit()  
}

# Read files
lines <- readLines(pop.info)
geno <- read.table(geno.file, header = T, sep = "\t")
geno.info <- read.table(geno.info.file, header = T, sep = "\t")

#separate genotype file into 3 lists
trio <- c()
geno.trio <- list()
for (i in 1:length(lines)){
  line.elements <- unlist(strsplit(lines[i],"\t"))
  trio <- c(trio, line.elements[1])
  line.elements <- line.elements[2:length(line.elements)] #contain only sample names
  geno.trio.curr <- geno[,colnames(geno) %in% line.elements]
  rownames(geno.trio.curr) <- geno[,1]
  geno.trio[[i]] <- geno.trio.curr
}
p1.cnt.anc <- rowSums(geno.trio[[1]] == "0", na.rm = T)
p1.cnt.alt <- rowSums(geno.trio[[1]] == "1", na.rm = T)
p1.cnt.total <- rowSums(!is.na(geno.trio[[1]]), na.rm = T)
p1.cnt.total.na <- rowSums(is.na(geno.trio[[1]]))
p2.cnt.anc <- rowSums(geno.trio[[2]] == "0", na.rm = T)
p2.cnt.alt <- rowSums(geno.trio[[2]] == "1", na.rm = T)
p2.cnt.total <- rowSums(!is.na(geno.trio[[2]]), na.rm = T)
p2.cnt.total.na <- rowSums(is.na(geno.trio[[2]]))

#prepare the allele site for introgression test
p1.anc.ratio <- p1.cnt.anc / p1.cnt.total
p1.alt.ratio <- p1.cnt.alt / p1.cnt.total
p1.na.ratio <- p1.cnt.total.na / (p1.cnt.total + p1.cnt.total.na)
p2.anc.ratio <- p2.cnt.anc / p2.cnt.total
p2.alt.ratio <- p2.cnt.alt / p2.cnt.total
p2.na.ratio <- p2.cnt.total.na / (p2.cnt.total + p2.cnt.total.na)
diff.anc.ratio <- p1.anc.ratio - p2.anc.ratio
diff.alt.ratio <- p1.alt.ratio - p2.alt.ratio
diff.ratio.table <- cbind(p1.anc.ratio, p1.alt.ratio, p1.na.ratio, p2.anc.ratio, p2.alt.ratio, p2.na.ratio, diff.anc.ratio, diff.alt.ratio)
#diff.ratio.table <- diff.ratio.table[complete.cases(diff.ratio.table),]
diff.ratio.table <- as.data.frame(diff.ratio.table)
#filter for the missing rate (<20%)
diff.ratio.table <- diff.ratio.table[diff.ratio.table$p1.na.ratio < 0.2 & diff.ratio.table$p2.na.ratio < 0.2,]
#filter for the difference >= 80%
diff.ratio.table <- diff.ratio.table[abs(diff.ratio.table$diff.anc.ratio) >= 0.8 & abs(diff.ratio.table$diff.alt.ratio) >= 0.8,]
#assign genotype to the table
P1_genotype <- ifelse(diff.ratio.table$p1.anc.ratio < 0.5, 1, 0)
P2_genotype <- ifelse(diff.ratio.table$p2.anc.ratio < 0.5, 1, 0)
diff.ratio.table <- cbind(P1_genotype, P2_genotype, diff.ratio.table)

#handle the test population
child.geno <- geno.trio[[3]][rownames(geno.trio[[3]]) %in% rownames(diff.ratio.table),]
child.geno.name <- rownames(child.geno)
child.pos <- c()
child.chr <- c()
for (i in 1:length(child.geno.name)){
  line.elements <- unlist(strsplit(child.geno.name[i], "_"))
  child.pos <- c(child.pos, line.elements[1])
  child.chr <- c(child.chr, line.elements[2])
}
child.geno <- cbind(Chr = child.chr, Pos = child.pos, P1_geno = diff.ratio.table$P1_genotype, P2_geno = diff.ratio.table$P2_genotype, child.geno)
child.geno[,2:ncol(child.geno)] <- child.geno[,2:ncol(child.geno)] %>% mutate_if(is.character, as.numeric)
#child.geno$average <- rowMeans(child.geno[7:ncol(child.geno)], na.rm = T) #the last column is the average of the test population

#handle individuals in the test population
pop.table <- c()
for (i in 5:ncol(child.geno)){
  indv.table <- child.geno[,c(1:4,i)]
  indv.table$matches <- indv.table[,ncol(indv.table)] == indv.table$P1_geno #determine the SNP identity, P1=TRUE, P2=FALSE
  indv.table$result <- ifelse(indv.table$matches, 0, 1) # assign 0 = P1; 1 = P2
  if (length(pop.table) == 0){
    pop.table <- indv.table[,c(1,2,7)]
  } else {
    pop.table <- cbind(pop.table, indv.table$result)
  }
  colnames(pop.table)[ncol(pop.table)] <- colnames(indv.table)[5]
  indv.unique.table <- indv.table %>%
    group_by(Chr,Pos) %>%
    summarise_at(c("result"), mean, na.rm = F)
  colnames(indv.unique.table)[3] <- colnames(indv.table)[5]
  indv.unique.table <- data.frame(indv.unique.table)
  if (i == 5){ #this is for introgression plot of all individuals
    pop.unique.table <- indv.unique.table
  } else {
    pop.unique.table <- cbind(indv.unique.table[,3])
  }
  chrs <- as.character(unlist(unique(indv.unique.table[,1]))) #detect how many chromosomes in the sample
  curr.unique.win <- data.frame(matrix(ncol = 3, nrow = 0))
  for (k in 1:length(chrs)){
    #get the length of the current chr
    chr.len <- as.numeric(geno.info[geno.info[,1] == chrs[k],2])
    indv.unique.table.chr <- indv.unique.table[indv.unique.table[,1] %in% chrs[k],]
    #remove the missing data in the individual
    indv.unique.table.chr <- indv.unique.table.chr[complete.cases(indv.unique.table.chr),]
    #duplicate the first and the last rows
    indv.unique.table.chr <- rbind(indv.unique.table.chr[1,],indv.unique.table.chr,indv.unique.table.chr[nrow(indv.unique.table.chr),])
    #change the position in the first and the last rows to match the actual chr length
    indv.unique.table.chr[1,2] <- 1
    indv.unique.table.chr[nrow(indv.unique.table.chr),2] <- chr.len
    #start to slide data into windows
    slide.number <- ceiling(nrow(indv.unique.table.chr) / step.size)
    start.win.pos <- c()
    n.win.start <- 1
    for (l in 1:slide.number){
      n.win.end = n.win.start + window.size - 1
      if (n.win.end > nrow(indv.unique.table.chr)){ #define the max border
        n.win.end = nrow(indv.unique.table.chr)
      }
      if (n.win.end - n.win.start < 10){ #if the last window has SNP number < 10, don't count it
        break
      }
      #focus on the current window
      slide.table <- indv.unique.table.chr[n.win.start:n.win.end,]
      slide.win.pos.avg <- mean(slide.table$Pos)
      slide.win.group.avg <- mean(slide.table[,ncol(slide.table)])
      slide.win.group.ci <- c()
      if (length(unique(slide.table[,ncol(slide.table)])) > 1){
        slide.win.group.ci <- t.test(slide.table[,ncol(slide.table)])$conf.int
      } else {
        slide.win.group.ci <- c(unique(slide.table[,ncol(slide.table)]),unique(slide.table[,ncol(slide.table)]))
      }
      win.row <- c(chrs[k], slide.win.pos.avg, slide.win.group.avg, slide.win.group.ci[2], slide.win.group.ci[1])
      curr.unique.win <- rbind(curr.unique.win, win.row)
      #focus on the current step window
      #message(paste0("debug: start ", n.step.start, " end ", n.step.end))
      colnames(curr.unique.win) <- c("Chr","Pos",colnames(indv.unique.table.chr)[ncol(indv.unique.table.chr)], "CI_95_lower", "CI_95_upper")
      n.win.start <- n.win.start + step.size
    }
  }
  curr.unique.win[,2:5] <- curr.unique.win[,2:5] %>% mutate_if(is.character, as.numeric)
  curr.unique.win[,3:5] <- 1 - curr.unique.win[,3:5] #ratio is actually opposite to the numbers
  curr.unique.win$CI_95_lower[curr.unique.win$CI_95_lower < 0] <- 0
  curr.unique.win$CI_95_upper[curr.unique.win$CI_95_upper > 1] <- 1
  if (length(chrs) == 1){
    save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[3], "_", chrs,"_introgression.rda"))
  } else {
    save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[3],"_introgression.rda"))
  }
}

#define 95% CI function
get_ci <- function(row) {
  if (length(unique(row)) == 1) {
    return(c(unique(row), unique(row))) # Return NA if the row is constant
  } else {
    return(t.test(row, conf.level = 0.95)$conf.int)
  }
}

#handle population table
pop.unique.table <- pop.table %>%
  group_by(Chr, Pos) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = 'keep')
pop.unique.table <- data.frame(pop.unique.table)
#average all the sample values
pop.unique.table$all <- rowMeans(pop.unique.table[,3:ncol(pop.unique.table)], na.rm = T)
colnames(pop.unique.table)[ncol(pop.unique.table)] <- trio[3]
curr.unique.win <- data.frame(matrix(ncol = 3, nrow = 0))
for (k in 1:length(chrs)){
  #get the length of the current chr
  chr.len <- as.numeric(geno.info[geno.info[,1] == chrs[k],2])
  pop.unique.table.chr <- pop.unique.table[pop.unique.table[,1] %in% chrs[k],]
  #duplicate the first and the last rows
  pop.unique.table.chr <- rbind(pop.unique.table.chr[1,],pop.unique.table.chr,pop.unique.table.chr[nrow(pop.unique.table.chr),])
  #change the position in the first and the last rows to match the actual chr length
  pop.unique.table.chr[1,2] <- 1
  pop.unique.table.chr[nrow(pop.unique.table.chr),2] <- chr.len
  slide.number <- ceiling(nrow(pop.unique.table.chr) / step.size)
  start.win.pos <- c()
  n.win.start <- 1
  for (l in 1:slide.number){
    n.win.end = n.win.start + window.size - 1
    if (n.win.end > nrow(pop.unique.table.chr)){ #define the max broader
      n.win.end = nrow(pop.unique.table.chr)
    }
    if (n.win.end - n.win.start < 10){ #if the last window has SNP number < 10, don't count it
      break
    }
    #focus on the current window
    slide.table <- pop.unique.table.chr[n.win.start:n.win.end,]
    slide.win.pos.avg <- mean(slide.table$Pos)
    slide.win.group.avg <- mean(slide.table[,ncol(slide.table)])
    if (length(unique(slide.table[,ncol(slide.table)])) > 1){
      slide.win.group.ci <- t.test(slide.table[,ncol(slide.table)])$conf.int
    } else {
      slide.win.group.ci <- c(unique(slide.table[,ncol(slide.table)]),unique(slide.table[,ncol(slide.table)]))
    }
    win.row <- c(chrs[k], slide.win.pos.avg, slide.win.group.avg, slide.win.group.ci[2], slide.win.group.ci[1])
    curr.unique.win <- rbind(curr.unique.win, win.row)
    #focus on the current step window
    #message(paste0("debug: start ", n.step.start, " end ", n.step.end))
    colnames(curr.unique.win) <- c("Chr","Pos",colnames(pop.unique.table.chr)[ncol(pop.unique.table.chr)], "CI_95_lower", "CI_95_upper")
    n.win.start <- n.win.start + step.size
  }
}
chrs <- as.character(unlist(unique(curr.unique.win[,1]))) #detect how many chromosomes in the sample
curr.unique.win[,2:5] <- curr.unique.win[,2:5] %>% mutate_if(is.character, as.numeric)
curr.unique.win[,3:5] <- 1 - curr.unique.win[,3:5] #ratio is actually opposite to the numbers
curr.unique.win$CI_95_lower[curr.unique.win$CI_95_lower < 0] <- 0
curr.unique.win$CI_95_upper[curr.unique.win$CI_95_upper > 1] <- 1
if (length(chrs) == 1){
  save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[3], "_", chrs,"_introgression.rda"))
} else {
  save(curr.unique.win, file = paste0(path,"/", colnames(curr.unique.win)[3],"_introgression.rda"))
}
