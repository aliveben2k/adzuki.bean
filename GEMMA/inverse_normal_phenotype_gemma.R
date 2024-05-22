#Usage: inverse_normal_phenotype.R phenotype_file
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
file <- as.character(args[1])
table <- read.table(file, sep = "\t", header = F)
for (i in 1:length(table)){
  uni.trait <- unique(na.omit(table[,i]))
  if (length(uni.trait) == 2){ #do not do inverse normal with binary data
    if (!grepl("0|1", uni.trait[1]) || !grepl("0|1", uni.trait[2])){ #convert binary data to 0/1 system
      table[,i] <- gsub(uni.trait[1], 99999999, table[,i])
      table[,i] <- gsub(uni.trait[2], 99999998, table[,i])
      table[,i] <- gsub(99999999, 0, table[,i])
      table[,i] <- gsub(99999998, 1, table[,i])
    }
  } else { #do inverse normal
    yranks <- rank(table[,i], na.last = "keep") #do not rank "na" value
    tempp <- (yranks-.5) / (length(yranks[!is.na(yranks)])) #do not count "na"
    table[,i] <- qnorm(tempp)
  }
}
outfile <- sub("pheno$", "inverse.pheno", file)
write.table(table, outfile, sep = "\t", quote = F, row.names = F, col.names = F)
