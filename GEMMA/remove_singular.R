args <- commandArgs(trailingOnly = TRUE)
pheno_col <- as.integer(args[1])
pheno_file <- as.character(args[2])
cov_file <- as.character(args[3])

pheno <- read.csv(pheno_file, sep = "\t", header = F)
cov <- read.csv(cov_file, sep = "\t", header = F)

pheno.exist <- which(!is.na(pheno[,pheno_col]))

all.uni.cov <- c()
for (i in 2:ncol(cov)){
  current.cov.values <- cov[pheno.exist, i]
  uni.cov <- unique(current.cov.values)
  if (length(uni.cov) == 1){
    all.uni.cov <- c(all.uni.cov, i)
  }
}

if (length(all.uni.cov) > 0){
  cov <- cov[,-all.uni.cov]
  if (ncol(cov) > 1){
    cov_file <- paste0(cov_file, '.tmp', pheno_col)
    write.table(cov, file = cov_file, sep = "\t", row.names = F, col.names = F, quote = F)
  }
  else {
    message("F")
  }
}
