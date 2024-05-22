library(data.table)
library(dplyr)
os <- Sys.info()[['sysname']]
#file.list.out <- c()
args <- commandArgs(trailingOnly = TRUE)
mu <- as.numeric(args[1])
gen <- as.numeric(args[2])
args[3] <- sub("/$|\\\\$", "", args[3])
dirs <- list.dirs(args[3], full.names = T, recursive = T)
path <- args[3]
files <- c()
pop_names <- c()
for (k in 1:length(dirs)){
  #get the file list
  files.tmp <- list.files(dirs[k], pattern="final\\.txt$", full.names=TRUE)
  files.tmp <- files.tmp[grepl("combined", files.tmp)]
  if (length(files.tmp) == 0){
    next
  }
  files <- c(files, files.tmp)
  #get the population name
  tmp <- unlist(strsplit(files.tmp[1],split="/|\\\\"))
  folder.name <- tmp[length(tmp)-1]
  count.underline <- lengths(regmatches(folder.name, gregexpr("_", folder.name)))
  pop.sep.underline <- ceiling(count.underline/2)
  tmp <- unlist(strsplit(folder.name,split="_"))
  if (length(tmp) == 2){
    pop.name.1 <- tmp[1]
    pop.name.2 <- tmp[2]
  } else {
    pop.name.1 <- paste0(tmp[1:pop.sep.underline],collapse = "_")
    pop.name.2 <- paste0(tmp[(pop.sep.underline+1):length(tmp)],collapse = "_")
  }
  pop_names <- c(pop_names, pop.name.1, pop.name.2)
  pop_names <- unique(pop_names)
}
if (length(files) == 0){
  quit("no")
}
f_order <- c()
for (l in 1:length(files)){ #get the file names without path
  tmp <- unlist(strsplit(files[l],split="/|\\\\"))
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
#do 2 by 2 populations averaging
for (n in 1:(length(pop_names)-1)){
  for (o in (n+1):length(pop_names)){
    current.files <- files[grepl(pop_names[n], files)]
    current.files <- current.files[grepl(pop_names[o], current.files)]
    #do calculation of each combination set
    ##combine all data in different runs together
    tmp_f <- c()
    for (j in 1:length(current.files)){
      possibleError <- tryCatch({
        tmp_table <- read.table(current.files[j], header = T, sep = "\t")}, 
        error = function(e) e)
      if(inherits(possibleError, "error")){ 
        message(paste0(current.files[j], ' is empty. Skip.'))
        next
      }
      left_time_boundary <- log10(tmp_table$left_time_boundary/mu*gen)
      center_time <- (log10(tmp_table$left_time_boundary/mu*gen)+log10(tmp_table$right_time_boundary/mu*gen))/2
      right_time_boundary <- log10(tmp_table$right_time_boundary/mu*gen)
      y <- (2*tmp_table$lambda_01)/(tmp_table$lambda_00+tmp_table$lambda_11)
      tmp_table <- as.data.frame(cbind(time_index = tmp_table$time_index, left_time_boundary, center_time, right_time_boundary, y))
      #message(tmp_table[2,])
      colnames(tmp_table) <- c("time_index","left_time_boundary", "center_time", "right_time_boundary", "y")
      tmp_table$y <- tmp_table$y/max(tmp_table$y, na.rm = T) #scale y value, so that maximum of the y value is 1
      tmp_f <- rbind(tmp_f, tmp_table)
    }
    info.data <- tmp_f %>% group_by(time_index) %>%
      summarise(mean.x.left = mean(left_time_boundary, na.rm = T),
                mean.x.center = mean(center_time, na.rm = T),
                mean.x.right = mean(right_time_boundary, na.rm = T), 
                mean.y = mean(y, na.rm = T),
                se.y = sd(y, na.rm = T)/sqrt(sum(!is.na(y) & y != 0)),
                n = sum(!is.na(y) & y != 0))
    #calculate 95% CI
    alpha <- 0.05
    degrees.freedom <- info.data$n[1] - 1
    t.score <- qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
    margin.error <- t.score * info.data$se.y
    upper.bound <- info.data$mean.y + margin.error
    lower.bound <- info.data$mean.y - margin.error
    info.data <- cbind(info.data,upper.bound,lower.bound)
    #path.out <- sub("/$","", path)
    #path.out <- paste0(path.out,"/",pop.names[i],".all.test.final.out")
    #write.table(info.data, file = path.out, quote = F, sep = "\t", row.names = F, col.names = T)  
    output.file <- unlist(strsplit(current.files[1],split="/|\\\\"))
    output.file <- paste(output.file[1:(length(output.file)-2)], sep = '/', collapse = '/')
    output.name <- paste0(pop_names[n], '_', pop_names[o], '.combined.all.final.out')
    output.file <- paste0(output.file, '/', output.name)
    if (grepl("windows", os, ignore.case=T)){
      output.file <- gsub("/", "\\\\", output.file)
    }
    write.table(info.data, file = output.file, quote = F, sep = "\t", row.names = F, col.names = T)  
    #file.list.out <- c(file.list.out, output.file)
  }
}
message("done")