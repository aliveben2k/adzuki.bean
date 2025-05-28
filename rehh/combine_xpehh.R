library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
window_size <- 10000
path <- sub("/$","", path)
files <- list.files(path, pattern="_xpehh.rda$", full.names=TRUE)
groups <- c()
for (i in 1:length(files)){
  strings <- unlist(strsplit(files[i], "\\."))
  strings[length(strings)-1] <- sub("_xpehh", "", strings[length(strings)-1])
  groups <- c(groups, strings[length(strings)-1])
}
groups <- unique(groups)
for (i in 1:length(groups)){
  group.files <- files[grepl(paste0(groups[i], "_xpehh"), files)]
  group.files <- sort(group.files)
  group.xpehh.windowed <- c()
  for (j in 1:length(group.files)){
    load(group.files[j])
    chr_length <- max(xp.ehh$POSITION)
    xp.ehh <- xp.ehh[complete.cases(xp.ehh),] #discard rows with a NA value
    xp.ehh$window <- cut(
      xp.ehh$POSITION,
      breaks = seq(0, chr_length, by = window_size),
      labels = seq(window_size, chr_length, by = window_size)
    )
    xp.ehh.windowed <- xp.ehh %>%
      group_by(window) %>%
      summarise(log10_p = mean(LOGPVALUE, na.rm = TRUE))
    #group.xpehh.windowed <- data.frame(group.xpehh.windowed)
    start.pos <- as.numeric(levels(xp.ehh.windowed$window))[xp.ehh.windowed$window] - window_size + 1
    xp.ehh.windowed$window <- as.numeric(levels(xp.ehh.windowed$window))[xp.ehh.windowed$window]
    group.name <- colnames(xp.ehh)[c(3,5)]
    group.name <- sub("FREQ_A.","", group.name)
    group.name.rep1 <- rep(group.name[1], nrow(xp.ehh.windowed))
    group.name.rep2 <- rep(group.name[2], nrow(xp.ehh.windowed))
    chr.name.rep <- rep(unique(xp.ehh$CHR), nrow(xp.ehh.windowed))
    xp.ehh.windowed <- cbind(pop1 = group.name.rep1, pop2 = group.name.rep2, CHR = chr.name.rep, start = start.pos, xp.ehh.windowed)
    colnames(xp.ehh.windowed)[c(5,6)] <- c("end","xpehh")
    group.xpehh.windowed <- rbind(group.xpehh.windowed, xp.ehh.windowed)
  }
  group.name <- paste0(group.name, collapse = "_")
  save(group.xpehh.windowed, file = paste0(path,"/", group.name, ".all.xpehh.windowed.rda"))
}