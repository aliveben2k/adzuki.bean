library(tidyverse)
library(splines)
library(ggplot2)
options(scipen = 999) #don't use scientific expression
args <- commandArgs(trailingOnly = TRUE)
recomb_file <- as.character(args[1]) #recomb file
path <- as.character(args[2]) # output folder path
no_sp <- as.integer(args[3]) # 0 = spline; 1 = no_spline
df_num <- as.integer(args[4]) # df of spline
prefix <- as.character(args[5]) # prefix of the contig name

ori <- read.table(recomb_file, header = T, sep = "\t")
chrs <- unique(ori[,1]) #get all unique chr names in the chr column
if (exists(prefix) || prefix != 0){ #filter for chr names that match the prefix
  chrs <- grep(paste0("^",prefix), chrs, value = T)
}
for (y in (1:length(chrs))){
  chr <- subset(ori, ori[,1]==chrs[y]) #subset the data by chr name
  chr[length(chr[,3]),3] <- chr[length(chr[,3])-1,3] #replace the last rate value with the previous value
  if (no_sp == 0){
    sspline <- smooth.spline(chr[,2], chr[,3], df=df_num) #make a spline
    tiff(paste0(path, '/', chrs[y], ".spline.imap.tiff"), units = "in", pointsize = 12, res = 300, bg = "white", compression = c("none"), width = 8, height = 5)
    plot(chr[,2], chr[,3], xlab="Position", ylab="Recombination Rate")
    lines(sspline, col="red")
    dev.off()
    out <- predict(sspline, x=chr[,2]) #grep spline value based on position
    table <- data.frame(matrix(unlist(out), ncol = length(out)))
    colnames(table) = c("pos","rate")
    cM <- rep(0, length(table$pos))
    table <- as_tibble(data.frame(table, cM))
    table[1,2] <- round(table[1,2],4)
    if (table[1,2] < 0){ #make sure rate value is bigger than 0 in the first row
      table[1,2] <- 0
    }
    for(x in 2:length(table$cM)){
      if (table[x,2] < 0){ #make sure rate value is bigger than 0 in the other row
        table[x,2] <- 0
      }
      table[x,2] <- round(table[x,2],4)
#      table[x,3] <- chr[x,4] #this code won't change cM value
      #the following code will change cM value based on new rate
      table[x,3] <- table[x-1,2]*(table[x,1]-table[x-1,1])/1000000+table[x-1,3] #re-calculate cM value
      table[x,3] <- round(table[x,3],4)
      if (table[x,3] < table[x-1,3]){ #make sure cM value is always bigger or equal to the [n-1] row value
       table[x,3] <- table[x-1,3]
      }
    }
    write.table(table, file=paste0(path, '/', chrs[y], ".spline.imap") , sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) #write table with impute5 format
    b_table <- as_tibble(data.frame(chr[,1],".",table[,3],chr[,2]))
    write.table(b_table, file=paste0(path, '/', chrs[y], ".spline.bmap") , sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) #write table with beagle5 format
  }
  if (no_sp == 1){
    table <- as_tibble(data.frame(chr[,2],chr[,3],chr[,4]))
    colnames(table) = c("pos","rate","cM")
    write.table(table, file=paste0(path, '/', chrs[y], ".imap") , sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) #write table with impute5 format
    b_table <- as_tibble(data.frame(chr[,1],".",chr[,4],chr[,2]))
    write.table(b_table, file=paste0(path, '/', chrs[y], ".bmap") , sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) #write table with beagle5 format
  }
}