options(warn=-1)
s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. May. 2023
Usage: Rscript rehh_calc.R -thap THAP_FILE -map MAP_FILE -chr CHR_NAME [-pos POSITION] [-pinfo INFO_FILE] [-popi POPULATION_NAMES] [-popi2 POPULATION_NAMES] [-syn] [-min MIN_BOUNDARY] [-max MAX_BOUNDARY] [-l EHH_LIMIT_VALUE] [-w WINDOW_SIZE] [-maf MIN_MAF_VALUE] [-gm GENO_MISSING_PREC] [-hm HAPLOTYPE_MISSING_PREC] [-m METHOD] [-h]\n
-thap: thap format file.
-map: map format file.
-chr: chromosome/contig name.
-pos: target position. (Numeric) (region eHH)
-pinfo: population information. (only used in ehhs and xpehh)
-popi: population(s) of interests. (only used in ehhs and xpehh)
-popi2: population(s) of interests. (the comparison pair of -popi, only used in xpehh)
-syn: indicating a synthetic F1 dataset. (only used in ehh and xpehh)
-min: left boundary of the culculation position (bp). Necessary when -pos is set.
-max: right boundary of the culculation position (bp). Necessary when -pos is set.
-l: threshold value for EHH(S) calculation. (0-1) Default: 0.01
-w: window size for calculation. (genome iHS) Default: 10000
-maf: cutoff value for MAF. Default: 0
-gm: missing percentage of the genotyped markers. Default: 0.15
-hm: missing percentage of the haplotypes on genotyped markers. Default: 0.1
-m: method to use. (ehh,ehhs,xpehh,ihs) Default: ehh
-h: help.\n\n")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}

thap_input <- c(); map_input <- c(); chr <- c(); w_size <- c(); rehh.method <- "ehh"
min.expand.range <- c(); max.expand.range <- c(); sp.mrk <- c(); popi <- c(); popi2 <- c()
lim.ehh <- c(); min.maf <- c(); mrk.missing <- c(); geno.missing <- c()
pinfo <- c(); syn <- 0;
for (i in 1:length(args)){
  if (args[i] == '-thap'){
    thap_input <- args[i+1]
    if (!file.exists(thap_input)){
      s.help()
      cat("-thap: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-map'){
    map_input <- args[i+1]
    if (!file.exists(map_input)){
      s.help()
      cat("-map: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-chr'){
    chr <- as.character(args[i+1])
  }
  if (args[i] == '-w'){
    w_size <- args[i+1]
    if (grepl("[^0-9]", w_size)){
      s.help()
      cat("-w: value is wrong.\n")
      quit()
    }
    w_size <- as.integer(w_size)
  }
  if (args[i] == '-min'){
    min.expand.range <- args[i+1]
    if (grepl("[^0-9]", min.expand.range)){
      s.help()
      cat("-min: value is wrong.\n")
      quit()
    }
    min.expand.range <- as.integer(min.expand.range)
  }
  if (args[i] == '-max'){
    max.expand.range <- args[i+1]
    if (grepl("[^0-9]", max.expand.range)){
      s.help()
      cat("-max: value is wrong.\n")
      quit()
    }
    max.expand.range <- as.integer(max.expand.range)
  }
  if (args[i] == '-pos'){
    sp.mrk <- args[i+1]
    if (grepl("[^0-9]", sp.mrk)){
      s.help()
      cat("-pos: value is wrong.\n")
      quit()
    }
    sp.mrk <- as.integer(sp.mrk)
  }
  if (args[i] == '-l'){
    lim.ehh <- args[i+1]
    if (grepl("[^0-9.]", lim.ehh)){
      s.help()
      cat("-l: value is wrong.\n")
      quit()
    }
    if (lim.ehh > 1 || lim.ehh < 0){
      s.help()
      cat("-l: value should be between 0 and 1.\n")
      quit()      
    }
    lim.ehh <- as.numeric(lim.ehh)
  }
  if (args[i] == '-maf'){
    min.maf <- args[i+1]
    if (grepl("[^0-9.]", min.maf)){
      s.help()
      cat("-maf: value is wrong.\n")
      quit()
    }
    if (min.maf > 1 || min.maf < 0){
      s.help()
      cat("-maf: value should be between 0 and 1.\n")
      quit()      
    }
    min.maf <- as.numeric(min.maf)
  }
  if (args[i] == '-hm'){ #missing percentage of the haplotype on genotyped markers
    mrk.missing <- args[i+1]
    if (grepl("[^0-9.]", mrk.missing)){
      s.help()
      cat("-hm: value is wrong.\n")
      quit()
    }
    if (mrk.missing > 1 || mrk.missing < 0){
      s.help()
      cat("-hm: value should be between 0 and 1.\n")
      quit()      
    }
    mrk.missing <- 100 - (mrk.missing * 100)
    mrk.missing <- as.integer(mrk.missing)
  }
  if (args[i] == '-gm'){ #missing percentage of the genotypes
    geno.missing <- args[i+1]
    if (grepl("[^0-9.]", geno.missing)){
      s.help()
      cat("-gm: value is wrong.\n")
      quit()
    }
    if (geno.missing > 1 || geno.missing < 0){
      s.help()
      cat("-gm: value should be between 0 and 1.\n")
      quit()      
    }
    geno.missing <- 100 - (geno.missing * 100)
    geno.missing <- as.integer(geno.missing)
  }
  if (args[i] == '-m'){
    rehh.method <- args[i+1]
    if (!grepl("\\behh\\b|\\behhs\\b|\\bxpehh\\b|\\bihs\\b", rehh.method)){
      s.help()
      cat("-m: wrong value. possible values: ehh,ehhs,xpehh,ihs\n")
      quit()
    }
  }
  if (args[i] == '-pinfo'){
    pinfo <- args[i+1]
    if (!file.exists(pinfo)){
      s.help()
      cat("-pinfo: file does not exist.\n")
      quit()
    }    
  }
  if (args[i] == '-popi'){
    popi <- args[i+1]
    if (popi == 'all'){
      popi <- c()
    }
  }
  if (args[i] == '-popi2'){
    popi2 <- args[i+1]
    if (popi2 == 'all'){
      popi2 <- c()
    }
  }
  if (args[i] == '-syn'){
    syn <- 1
  }
  if (args[i] == '-h'){
    s.help()
    quit()
  }
}

if(!require("rehh")) install.packages("rehh")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("tidyverse")) install.packages("tidyverse")
library(rehh)
library(ggplot2)
library(tidyverse)
col.brewer.theme <- c("#4b8bcb","#ed3325","#7758a5","#255271", "#f7931e", "#921a1d", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")

if (length(thap_input) < 1){
  s.help()
  cat("-thap: argument must be provided.\n")
  quit()
}
if (length(map_input) < 1){
  s.help()
  cat("-map: argument must be provided.\n")
  quit()
}
if (length(lim.ehh) == 0){
  lim.ehh <- 0.01
}
if (length(min.maf) == 0){
  min.maf <- NA
}
x_min <- c(); x_max <- c()
if (length(min.expand.range) != 0){
  x_min = min.expand.range
  x_max = max.expand.range
}
if (length(geno.missing) == 0){
  geno.missing <- 85
}
if (length(mrk.missing) == 0){
  mrk.missing <- 90
}

pname <- c()
if (length(popi) != 0){
  pname <- gsub(',','_',popi)
} else {
  pname <- 'all'
}
if (length(popi2) != 0){
  pname2 <- gsub(',','_',popi2)
  pname <- paste0(pname,'_', pname2)
}

out <- sub(".thap", ".", thap_input)
if (file.exists(paste0(out, "td_", lim.ehh, ".", pname, ".haplohh.rda"))){
  load(paste0(out, "td_", lim.ehh, ".", pname, ".haplohh.rda"))
} else {
  gc()
  #transform vcf file and map file into haplohh format
  hap <- data2haplohh(hap_file = thap_input, map_file = map_input, haplotype.in.columns = TRUE, chr.name = chr, remove_multiple_markers = FALSE, allele_coding = "map", min_maf = NA, min_perc_geno.mrk = mrk.missing, min_perc_geno.hap = NA) #min_perc_geno.hap = geno.missing
  gc()
  #save data
  save(hap, file = paste0(out, "td_", lim.ehh, ".", pname, ".haplohh.rda"))
}
if (length(sp.mrk) != 0 && rehh.method == "ehh"){
  sp.mrk <- paste0(chr, "_", sp.mrk)
  if (length(min.expand.range) == 0 || length(max.expand.range) == 0){
    s.help()
    cat("-min/-max: argument must be provided when using -pos.\n")
    quit()
  }
  message("Calculating eHH of the region...")
  #calculate ehh of specific marker region
  hap.filter <- subset(hap, min_perc_geno.hap = geno.missing, min_maf = min.maf)
  sp.ehh.sp.mrk <- calc_ehh(hap.filter, sp.mrk, limhaplo = 2, limhomohaplo = 2, limehh = lim.ehh, phased = TRUE, polarized = TRUE, include_zero_values = TRUE)
  center_pos <- hap.filter@positions[[sp.mrk]]/1000000
  message("Calculating bi-furcation of the region...")
  #calculate bi-furcation of specific marker region
  furc <- calc_furcation(hap.filter, sp.mrk, allele = NA, limhaplo = 2, phased = TRUE, polarized = TRUE)
  #plot EHH for specific region
  plot.data <- sp.ehh.sp.mrk$ehh %>%
    mutate(pos = POSITION/1000000)
  ehh.all <- c()
  for (j in 2:(ncol(plot.data)-1)){
    curr.data <- plot.data[,c(ncol(plot.data),j)]
    type.value <- c()
    if (grepl("EHH_A",colnames(plot.data)[j])){
      type.value <- sub("EHH_A", "Ancestral", colnames(plot.data)[j])
    } else {
      type.value <- sub("EHH_D", "Derived", colnames(plot.data)[j])
    }
    curr.data <- cbind(curr.data, rep(type.value, nrow(plot.data)))
    colnames(curr.data) <- c("pos","EHH","type")
    ehh.all <- rbind(ehh.all, curr.data)
  }
  #save data
  save(sp.ehh.sp.mrk, ehh.all, furc, center_pos, sp.mrk, file = paste0(out, "td_", lim.ehh, ".", sp.mrk, ".ehh.rda"))
  #plot
  ehh.plot <- ggplot(ehh.all, aes(pos,EHH,color=type)) +
    geom_vline(xintercept = center_pos, linetype="dotted", color = "black", linewidth = 0.2) +
    geom_line(linewidth = 0.5) +
    theme_minimal() + 
    theme(aspect.ratio = 1,
          #axis.line = element_line(colour = "black", size = 0.2),
          axis.ticks = element_line(colour = "black", size = 0.2),
          text = element_text(size = 7.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 7.5),
          legend.spacing.y = unit(-0.2, 'cm'),
          axis.text = element_text(color = "black", size = 7.5),
          axis.title.x = element_text(color = "black", size = 7.5),
          axis.title.y = element_text(color = "black", size = 7.5),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          panel.grid = element_blank()
    ) +
    #ylim(0,1) +
    scale_x_continuous(name=paste0(hap.filter@chr.name, ' (Mb)'), limits = c(x_min/1000000,x_max/1000000)) +
    scale_y_continuous(name="EHH", limits =c(0,1)) +
    scale_color_manual(values = col.brewer.theme)
    #labs(x=paste0(hap.filter@chr.name, ' (Mb)'), y="EHH")
  tiff_out <- c()
  if (length(min.maf) > 0){
    tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", "maf_", min.maf, ".ehh.tiff")
  } else {
    tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".ehh.tiff")
  }
  tiff(tiff_out, units = "cm", res = 600, width = 8, height = 4)
  print(ehh.plot)
  dev.off()
  #plot bi-furcation for specific region
  for (i in 1:length(sp.ehh.sp.mrk$freq)){
    if (sp.ehh.sp.mrk$freq[i] == 0){
      next
    }
    out_name <- names(sp.ehh.sp.mrk$freq)[i]
    out_name <- gsub("FREQ","furcation", out_name)
    if (length(min.maf) > 0){
      tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", "maf_", min.maf, ".", out_name,".tiff")
    } else {
      tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", out_name,".tiff")
    }
    tiff(paste0(out, sp.mrk, ".", out_name,".tiff"), units = "cm", res = 600, width = 8, height = 4)
    plot(furc, allele = i-1, col = col.brewer.theme, mrk.col = "gray", lwd = 0.1, hap.names = NULL, cex.lab = 1, family.lab = "sans", offset.lab = 0.5, legend = NA, legend.xy.coords = "none")
    dev.off()
  }
} else if (rehh.method == "ihs") {
  if (length(w_size) == 0){
    w_size <- 10000
  }
  message("Calculating iHS of the genome...")
  #scan ihh of the genome
  scan.ihh <- scan_hh(hap, limhaplo = 2, limhomohaplo = 2, limehh = lim.ehh, limehhs = lim.ehh, phased = TRUE, polarized = TRUE)
  #convert ihh to ihs
  scan.ihs <- ihh2ihs(scan.ihh, freqbin = 0, min_maf = min.maf)
  #calculate ihs by window
  scan.ihs.window <- calc_candidate_regions(scan.ihs, window_size = w_size, pval = TRUE, threshold = 0)
  #save data without sp.ehh
  save(scan.ihh, scan.ihs, scan.ihs.window, file = paste0(out, lim.ehh, "_", chr, ".rda"))
  # #draw distribution plot
  # if (!file.exists(paste0(out, "td_", lim.ehh, ".distrib_plot.tiff"))){
  #   tiff(paste0(out, "td_", lim.ehh,".distrib_plot.tiff"), units = "cm", res = 600, width = 4, height = 4)
  #   distribplot(scan.ihs[["ihs"]][["IHS"]], lty = 1, lwd = 1, col = c("#4b8bcb", "#ed3325"), qqplot = FALSE)
  #   dev.off()
  # }
  # #manhattanplot
  # if (!file.exists(paste0(out, "td_", lim.ehh, ".manhat_plot.tiff"))){
  #   tiff(paste0(out, "td_", lim.ehh, ".manhat_plot.tiff"), units = "cm", width = 9, height = 4)
  #   manhattanplot(scan.ihs, pval = FALSE, threshold = c(-3, 3))
  #   dev.off()
  # }
} else if (length(sp.mrk) != 0 && rehh.method == "ehhs") {
  #col.brewer.theme <- c("#4b8bcb","gold","#ed3325","#255271", "#f7931e", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
  col.brewer.theme <- c("#4b8bcb","gold","#ed3325","#255271", "#c6b1d4", "#921a1d", "#7758a5", "#f9bbb9", "#c6b1d4", "#ed7f6d", "#90a7b7")
  sp.mrk <- paste0(chr, "_", sp.mrk)
  ehhs.pop <- c()
  if (length(pinfo) != 0){
    pinfo <- read.table(pinfo, header = F)
    pinfo.col <- c()
    if (ncol(pinfo) == 1){
      pinfo <- cbind(pinfo, rep('all', nrow(pinfo)))
      pinfo.col <- 2
    } else if (syn == 0){
      pinfo.col <- 2
    } else {
      pinfo.col <- 4
    }
    colnames(pinfo)[pinfo.col] <- "pop"
    if (length(popi) != 0){ #setup population if -popi exist
      ehhs.pop <- unlist(strsplit(popi, ','))
    } else { #calculate all populations one by one
      ehhs.pop <- unique(pinfo[,pinfo.col])
    }
  } else { #no sample list provided
    ehhs.pop <- "all"
  }
  ehhs.all <- c()
  for (i in 1:length(ehhs.pop)){
    if (ehhs.pop[i] != "all"){ #subset by population if pinfo exist
      curr.pinfo <- pinfo[pinfo$pop==ehhs.pop[i],]
      idx <- c(as.integer(rownames(curr.pinfo))*2 - 1, as.integer(rownames(curr.pinfo))*2)
      idx <- sort(idx)
      hap.subset <- subset(hap, select.hap = idx, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
    } else { #no pinfo, use all haps
      hap.subset <- subset(hap, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
    }
    center_pos <- hap.subset@positions[[sp.mrk]]/1000000
    curr.ehhs <- calc_ehhs(hap.subset, sp.mrk, limhaplo = 2, limhomohaplo = 2, limehhs = lim.ehh, include_zero_values = FALSE, include_nhaplo = TRUE, phased = TRUE)
    curr.data <- curr.ehhs$ehhs %>%
      mutate(pos = POSITION/1000000) %>%
      mutate(type = ehhs.pop[i])
    ehhs.all <- rbind(ehhs.all, curr.data)
  }
  save(ehhs.all, center_pos, ehhs.pop, pinfo, sp.mrk, file = paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", pname, ".ehhs.rda"))
  ehhs.all$type <- factor(ehhs.all$type, levels = c("JP_LR","CN_LRN","CN_LRS","JP_WL","SOUTH_LR"))
  ehhs.plot <- ggplot(ehhs.all, aes(pos,NEHHS,color=type)) +
    geom_vline(xintercept = center_pos, linetype="dotted", color = "black", linewidth = 0.2) +
    geom_line(linewidth = 0.5) +
    theme_minimal() +
    theme(aspect.ratio = 1,
          #axis.line = element_line(colour = "black", size = 0.2),
          axis.ticks = element_line(colour = "black", size = 0.2),
          text = element_text(size = 7.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 7.5),
          legend.spacing.y = unit(-0.2, 'cm'),
          axis.text = element_text(color = "black", size = 7.5),
          axis.title.x = element_text(color = "black", size = 7.5),
          axis.title.y = element_text(color = "black", size = 7.5),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          panel.grid = element_blank()
    ) +
    scale_x_continuous(name=paste0(hap.subset@chr.name, ' (Mb)'), breaks = seq(round(x_min/1000000, digits = 2),round(x_max/1000000, digits = 2),0.02)) +
    #scale_x_continuous(name=paste0(hap.subset@chr.name, ' (Mb)')) +
    scale_y_continuous(name="EHHS", limits =c(0,1)) +
    #ylim(0,1) +
    scale_color_manual(values = col.brewer.theme)
    #labs(x=paste0(hap.subset@chr.name, ' (Mb)'), y="EHHS")
  tiff_out <- c()
  if (length(min.maf) > 0){
    tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", "maf_", min.maf, ".", pname,".ehhs.tiff")
  } else {
    tiff_out <- paste0(out, "td_", lim.ehh, ".", sp.mrk, ".", pname,".ehhs.tiff")
  }
  tiff(tiff_out, units = "cm", res = 600, width = 8, height = 4)
  print(ehhs.plot)
  dev.off()    
} else if (rehh.method == "xpehh") {
  if (length(w_size) == 0){
    w_size <- 10000
  }
  message("Calculating XP-EHH of the genome...")
  ehhs.pop <- c()
  ehhs.pop2 <- c()
  if (length(pinfo) != 0){
    pinfo <- read.table(pinfo, header = F)
    pinfo.col <- c()
    if (ncol(pinfo) == 1){
      pinfo <- cbind(pinfo, rep('all', nrow(pinfo)))
      pinfo.col <- 2
    } else if (syn == 0){
      pinfo.col <- 2
    } else {
      pinfo.col <- 4
    }
    colnames(pinfo)[pinfo.col] <- "pop"
    if (length(popi) != 0){ #setup population if -popi exist
      ehhs.pop <- unlist(strsplit(popi, ','))
      ehhs.pop2 <- unlist(strsplit(popi2, ','))
    } else { #calculate all populations one by one
      ehhs.pop <- unique(pinfo[,pinfo.col])
      ehhs.pop2 <- unique(pinfo[,pinfo.col])
    }
  } else { #no sample list provided
    ehhs.pop <- "all"
    ehhs.pop2 <- "all"
  }
  hap.subset.all <- c()
  hap.subset2.all <- c()
  if (ehhs.pop[1] != "all"){ #subset by population if pinfo exist
    curr.pinfo <- pinfo[pinfo$pop %in% ehhs.pop,]
    idx <- c(as.integer(rownames(curr.pinfo))*2 - 1, as.integer(rownames(curr.pinfo))*2)
    idx <- sort(idx)
    hap.subset <- subset(hap, select.hap = idx, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
    curr.pinfo2 <- pinfo[pinfo$pop %in% ehhs.pop2,]
    idx2 <- c(as.integer(rownames(curr.pinfo2))*2 - 1, as.integer(rownames(curr.pinfo2))*2)
    idx2 <- sort(idx2)
    hap.subset2 <- subset(hap, select.hap = idx2, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
  } else { #no pinfo, use all haps
    hap.subset <- subset(hap, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
    hap.subset2 <- subset(hap, min_perc_geno.hap = geno.missing, min_maf = min.maf, min_perc_geno.mrk = mrk.missing)
  }
  #hap.subset.all <- rbind(hap.subset.all, hap.subset)
  #hap.subset.all2 <- rbind(hap.subset.all2, hap.subset2)
  
  #scan iES of the genome
  if (!file.exists(paste0(out, pname, "_", chr,"_scan_hh_xpehh.rda"))){
    scan.ies <- scan_hh(hap.subset, limhaplo = 2, limhomohaplo = 2, limehh = lim.ehh, limehhs = lim.ehh, phased = TRUE, polarized = TRUE)
    scan.ies2 <- scan_hh(hap.subset2, limhaplo = 2, limhomohaplo = 2, limehh = lim.ehh, limehhs = lim.ehh, phased = TRUE, polarized = TRUE)
    save(scan.ies, scan.ies2, lim.ehh, min.maf, w_size, ehhs.pop, ehhs.pop2, file = paste0(out, pname, "_", chr,"_scan_hh_xpehh.rda"))
  } else {
    load(paste0(out, pname, "_", chr,"_scan_hh_xpehh.rda"))
  }
  if (!file.exists(paste0(out, pname, "_", chr,"_xpehh.rda"))){
    #calculate XP-EHH
    xp.ehh <- ies2xpehh(scan.ies, scan.ies2, ehhs.pop, ehhs.pop2)
    #save data without sp.ehh
    save(xp.ehh, file = paste0(out, pname,"_xpehh.rda"))
  }
}


