# Rscript this_code.R input_file output_file
# Calculate pairwise genomic distance within each window
# This code directly accepts output from "vcf2table_missingNA.pl": first 5 columns are SNP info: name, chr, position, ref, alt
# Missing data in site calls are not correctly listed as NA
library(vroom)

args <- commandArgs(trailingOnly = TRUE)
input.file <- as.character(args[1])
output.file <- as.character(args[2])
message("reading txt file...")
if (grepl("tmp.txt",input.file,fixed=T)) {
  #data <- read.table(file=input.file,as.is=T,header=T,sep="\t")
  data <- tryCatch({
    if (file.size(input.file) > 0){
      vroom(input.file, delim = "\t", col_names = T)
    }
  }, error = function(err) {
    message(paste0(input.file," is empty.", sep = ""))
  })
} else {
  stop("Your input SNP matrix file should be in .txt format~!\n")
}
# data.filtered <- c()
# for (i in 1:nrow(data)){
#  test <- na.omit(unique(unlist(data[i,6:ncol(data)])))
#  if (length(test) > 1){
#    data.filtered <- rbind(data.filtered, data[i,])
#  }
# }
# data <- c()
# data <- data.filtered
message("file is read.")

sample.name <- names(data)[6:ncol(data)]

data.use <- as.matrix(data[,6:ncol(data)])
message("calculating...")
avail.list <- list()
diff.list <- list()
start_time <- Sys.time()
for (i in 1:(ncol(data.use)-1)) {
  #message(i)
  available.vec <- c()
  difference.vec <- c()
  for (j in (i+1):ncol(data.use)) {
    diff.vec.abs <- abs(as.numeric(data.use[,i]) - as.numeric(data.use[,j]))
    available <- sum(!is.na(diff.vec.abs))
    difference <- sum(diff.vec.abs,na.rm=T)
    available.vec <- c(available.vec,available)
    difference.vec <- c(difference.vec,difference)
  }
  avail.list[[i]] <- available.vec
  diff.list[[i]] <- difference.vec
}
available.vec <- unlist(avail.list)
difference.vec <- unlist(diff.list)
end_time <- Sys.time()
message(paste(i, end_time-start_time, sep = "\t"))
message("generating available matrix...")
available.matrix <- matrix(0,ncol=ncol(data.use),nrow=ncol(data.use))
available.matrix[lower.tri(available.matrix)] <- available.vec
t.available.matrix <- t(available.matrix)
available.matrix[upper.tri(available.matrix)] <- t.available.matrix[upper.tri(t.available.matrix)]

message("generating difference matrix...")
difference.matrix <- matrix(0,ncol=ncol(data.use),nrow=ncol(data.use))
difference.matrix[lower.tri(difference.matrix)] <- difference.vec
t.difference.matrix <- t(difference.matrix)
difference.matrix[upper.tri(difference.matrix)] <- t.difference.matrix[upper.tri(t.difference.matrix)]

#distance.matrix <- difference.matrix / available.matrix
#diag(distance.matrix) <- 0
#output.max <- paste(output.file,".mtx",sep="")
#row.names(distance.matrix) <- sample.name
#write.table(distance.matrix, file = output.max ,sep = "\t", quote = FALSE)

save(sample.name,available.matrix,difference.matrix,file=paste(output.file,".rda",sep=""))
message("done.")

