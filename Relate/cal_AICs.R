options(stringsAsFactors = FALSE)
files <- c(); seris <- c()
file.loci <- c()
path <- "."
files <- list.files(path, pattern="txt$", full.names=TRUE)
for (i in 1:length(files)){
  path.tmp <- unlist(strsplit(files[i],"/|\\\\"))
  seris <- path.tmp[length(path.tmp)]
  file.locus <- unlist(strsplit(files[i], "_CLUES_"))[2]
  file.locus <- sub("_inference\\.txt$","", file.locus) #no path
  file.loci <- c(file.loci, file.locus)
}
file.loci <- unique(file.loci)

all.inferences <- list()
for (i in 1:length(file.loci)){
  table.inference <- c()
  files.curr <- files[grepl(file.loci[i], files)]
  for (j in 1:length(files.curr)){
    tmp <- read.table(files.curr[j], header = T)
    tmp <- cbind(locus = file.loci[i], tmp)
    if (length(table.inference) == 0){
      table.inference <- tmp
    } else {
      table.inference <- rbind(table.inference, tmp)
    }
  }
  all.inferences[[i]] <- table.inference
}
AIC.table <- data.frame(locus = character(), AIC = numeric())
for (i in 1:length(all.inferences)){
  AICs.sum <- 0
  k.num <- (ncol(all.inferences[[i]])-3)/3
  for (j in 1:nrow(all.inferences[[i]])){
    #AIC values
    AIC.curr <- 2*k.num-2*all.inferences[[i]][j,2]
    AICs.sum <- AICs.sum+AIC.curr
    s.values <- list() #how many epochs, and how many s columns
    for (k in 1:k.num){
      is_present <- k >= 1 && k <= length(s.values)
      if (!is_present){
        s.values[[k]] <- all.inferences[[i]][j,(3+k*3)]
      } else {
        s.values[[k]] <- s.values[[k]]+all.inferences[[i]][j,(3+k*3)]
      }
    }
  }
  s.values.avg <- c()
  for (j in 1:k.num){
    s.values.avg.curr <- s.values[[j]]/nrow(all.inferences[[i]])
    s.values.avg <- c(s.values.avg, s.values.avg.curr)
  }
  AIC.avg <- AICs.sum/nrow(all.inferences[[i]])
  locus <- all.inferences[[i]][1,1]
  line <- c(locus, AIC.avg, s.values.avg)
  AIC.table <- rbind(AIC.table, line)
}
s.names <- seq(1,ncol(AIC.table)-2, by = 1)
s.names <- paste0("SelectionMLE", s.names)
colnames(AIC.table) <- c("locus","AIC", s.names)
write.table(AIC.table, file = "AIC_average.txt", quote = F, col.names = T, row.names = F, sep = "\t")