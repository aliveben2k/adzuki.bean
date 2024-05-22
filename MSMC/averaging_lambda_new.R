library(data.table)
library(dplyr)
file.list.out <- c()
args <- commandArgs(trailingOnly = TRUE)
dirs <- list.dirs(args[1], full.names = T, recursive = T)
path <- args[1]
files <- c()
for (k in 1:length(dirs)){
  files.tmp <- list.files(dirs[k], pattern="final\\.txt$", full.names=TRUE)
  files.tmp <- files.tmp[!grepl("combined", files.tmp)]
  tmp <- unlist(strsplit(files.tmp[1],split="/"))
  folder.name <- tmp[length(tmp)-1]
  files <- c(files, files.tmp[!grepl(paste0(folder.name, "_"), files.tmp)])
}
if (length(files) == 0){
  quit("no")
}
f_order <- c()
for (l in 1:length(files)){ #get the file names without path
  tmp <- unlist(strsplit(files[l],split="/"))
  f_order <- c(f_order, tmp[length(tmp)])
}
f_order.name <- f_order[order(f_order)] #reorder file names
f_order.idx <- order(f_order) #get the reordered index
files <- files[f_order.idx] #sort original file list with reordered index
name.table <- as.data.table(f_order.name)[, list(list(.I)), by = f_order.name] #make a name table with index (repeated files grouped together)
idx.unique <- c()
for (m in 1:nrow(name.table)){ #get the unique file index
  idx.unique <- c(idx.unique, unlist(name.table$V1[m])[1])
}
files <- files[idx.unique] #subset the original files with the unique file index
#start to do combination
current_pop <- c()
pop.data <- list()
pop.names <- c()
for (j in 1:length(files)){
  tmp.path <- unlist(strsplit(files[j],split="/"))
  tmp.name <- unlist(strsplit(tmp.path[length(tmp.path)],split="\\."))
  current_pop.tmp <- tmp.name[1]
  current_pop.tmp <- sub("_msmc2","",current_pop.tmp) #get the population name
  ##Va use
  if (lengths(strsplit(current_pop.tmp, '_')) > 2){
    next
  }
  ##Va use
  if (length(current_pop) == 0 || current_pop != current_pop.tmp){
    current_pop <- current_pop.tmp
    tmp_table <- read.table(files[j], header = T, sep = "\t")
    pop.data[[current_pop]] <- as.data.frame(tmp_table)
    pop.names <- c(pop.names, current_pop)
  } else {
    tmp_table <- read.table(files[j], header = T, sep = "\t")
    pop.data[[current_pop]] <- rbind(pop.data[[current_pop]], as.data.frame(tmp_table))
  }
}
for (i in 1:length(pop.names)){
  info.data <- pop.data[[pop.names[i]]] %>% group_by(time_index) %>% 
    summarise(mean.left = mean(left_time_boundary, na.rm = T), 
              mean.right = mean(right_time_boundary, na.rm = T), 
              mean.lambda = mean(lambda, na.rm = T), 
              se.lambda = sd(lambda, na.rm = T)/sqrt(sum(!is.na(lambda) & lambda != 0)),
              n = sum(!is.na(lambda) & lambda != 0))
  #calculate 95% CI
  alpha <- 0.05
  degrees.freedom <- info.data$n[1] - 1
  t.score <- qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  margin.error <- t.score * info.data$se.lambda
  upper.bound <- info.data$mean.lambda + margin.error
  lower.bound <- info.data$mean.lambda - margin.error
  info.data <- cbind(info.data,upper.bound,lower.bound)
  path.out <- sub("/$","", path)
  path.out <- paste0(path.out,"/",pop.names[i],".all.test.final.out")
  write.table(info.data, file = path.out, quote = F, sep = "\t", row.names = F, col.names = T)
}