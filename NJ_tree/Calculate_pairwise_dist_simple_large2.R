# Rscript this_code.R input_file output_file
# Calculate pairwise genomic distance within each window
# This code directly accepts output from "vcf2table_missingNA.pl": first 5 columns are SNP info: name, chr, position, ref, alt
# Missing data in site calls are not correctly listed as NA

args <- commandArgs(trailingOnly = TRUE)
path <- as.character(args[1])
output.file <- as.character(args[2])
input.files <- list.files(path, pattern='tmp.rda', full.names=TRUE)

available.matrix.all <- c()
difference.matrix.all <- c()
sample.name.all <- c()
for (i in 1:length(input.files)){
  load(input.files[i])
  if (i == 1){
    sample.name.all <- sample.name
    message(length(sample.name.all))
    available.matrix.all <- as.matrix(available.matrix)
    difference.matrix.all <- as.matrix(difference.matrix)
  } else {
    available.matrix.all = as.matrix(available.matrix.all) + as.matrix(available.matrix)
    difference.matrix.all = as.matrix(difference.matrix.all) + as.matrix(difference.matrix)
  }
  message(nrow(available.matrix.all))
}
distance.matrix <- as.matrix(difference.matrix.all / available.matrix.all)
diag(distance.matrix) <- 0
message(nrow(distance.matrix))
row.names(distance.matrix) <- sample.name.all
output.max <- paste(output.file,".mtx",sep="")
write.table(distance.matrix, file = output.max ,sep = "\t", quote = FALSE)

save(sample.name.all,available.matrix.all,difference.matrix.all,distance.matrix,file=paste(output.file,".rda",sep=""))

