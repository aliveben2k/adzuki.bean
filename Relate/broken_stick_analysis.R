# Load necessary libraries
library(segmented)
library(ggplot2)

file.open <- c()
path <- c()
npsi <- c()
for (i in 1:length(args)){
  if (args[i] == '-f'){
    file.open <- as.character(args[i+1]) #final_all_loci.broken.stick.rda
    path.tmp <- unlist(strsplit(file.open,"/|\\\\"))
    path <- paste0(path.tmp[1:(length(path.tmp)-1)], collapse = "/")
  }
  if (args[i] == '-n'){
    npsi <- as.numeric(args[i+1])
  }
}

# Load your data (replace 'your_data.csv' with your file)
load(file.open) #y.max.broken (loci list with repeat values), plot.data (final median table)
output.txt <- data.frame(locus = as.character(), npsi = as.numeric(), breakpoint = as.character(), breakpoint.min = as.character(), breakpoint.max = as.character())
for (i in 1:length(y.max.broken)){
  data <- y.max.broken[[i]] # get individual locus data
  # Get the current locus name
  locus.name <- unique(data$locus)
  data.median <- plot.data[plot.data$label == locus.name,]

  ## broken stick with median data
  time <- data.median$Time
  frequency <- data.median$Frequency
  lm_fit.median <- lm(frequency ~ time)
  for (npsi.num in 1:npsi){
    repeat.predicted_values <- list()
    segmented_fit.median <- segmented(lm_fit.median, seg.Z = ~ time, npsi = npsi.num)
    breakpoint.median <- segmented_fit.median$psi[,2]
    #breakpoint.min <- breakpoint.median
    #breakpoint.max <- breakpoint.median
    predicted_values.median <- predict(segmented_fit.median)
    ##
    # Time
    time <- data$time  # X-axis (e.g., years or generations)
    breakpoint.all <- c()
    for (j in 3:ncol(data)){ # Each repeat in a locus
      frequency <- data[,j]  # Y-axis (allele frequencies)
      # Fit a linear model (this is the base model)
      lm_fit <- lm(frequency ~ time)
      # Fit the segmented regression model with 1 breakpoint (npsi = 1)
      segmented_fit <- segmented(lm_fit, seg.Z = ~ time, npsi = npsi.num)
      # Get the estimated breakpoint
      breakpoint <- t(data.frame(segmented_fit$psi[,2]))
      if (j == 3){
        breakpoint.all <- breakpoint
      } else {
        
        breakpoint.all <- rbind(breakpoint.all, breakpoint)
      }
      predicted_values <- as.data.frame(predict(segmented_fit))
      predicted_values <- cbind(time = data$time, values = predicted_values)
      colnames(predicted_values) <- c("time", "value")
      repeat.predicted_values[[j]] <- predicted_values
      # Perform an ANOVA to compare the linear model with the segmented model
      # anova(lm_fit, segmented_fit)
      # Extract the slopes (selection coefficients) from the segmented model
      # slopes <- slope(segmented_fit)
    }
    all.remove.idx <- c()
    breakpoint.min <- c()
    breakpoint.max <- c()
    for (m in 1:ncol(breakpoint.all)){
      idx.top3 <- order(breakpoint.all[,m], decreasing = F)[1:3]
      idx.bottom3 <- order(breakpoint.all[,m], decreasing = T)[1:3]
      all.remove.idx <- c(all.remove.idx, idx.top3, idx.bottom3)
      breakpoint.min <- c(breakpoint.min, sort(breakpoint.all[,m], decreasing = F)[4])
      breakpoint.max <- c(breakpoint.max, sort(breakpoint.all[,m], decreasing = T)[4])
    }
    all.remove.idx <- unique(all.remove.idx)
    all.remove.idx <- sort(all.remove.idx, decreasing = T)
    all.remove.idx <- all.remove.idx + 2
    repeat.predicted_values <- repeat.predicted_values[-all.remove.idx]
    colors.palette <- c("deepskyblue3","darkgoldenrod3","darkorchid3")
    ##adzuki use
    color.plot <- c()
    adzuki.loci <- c("CN_LRNCN_LRSJP_LRSOUTH_LR_chr1_13815643-13815643", "CN_LRNCN_LRSJP_LRSOUTH_LR_chr4_4148520-4148520", "CN_LRNCN_LRSJP_LRSOUTH_LR_chr7_2492619-2492619")
    for (k in 1:length(adzuki.loci)){
      if (unique(locus.name) == adzuki.loci[k]){
        color.plot <- colors.palette[k]
      }
    }
    ##
    comp <- ggplot(data.median)
    for (l in 1:length(breakpoint.median)){
      comp <- comp +
        geom_rect(xmin = breakpoint.min[l], xmax = breakpoint.max[l], ymin = 0, ymax = 1, fill = "khaki1", alpha = 0.6)
    }
    for (k in seq_along(repeat.predicted_values)){
      if (length(repeat.predicted_values[[k]]) == 0){
        next
      }
      comp <- comp +
        geom_line(data = repeat.predicted_values[[k]], mapping = aes(time, value), color = "grey80", linewidth = 0.1, alpha = 0.8)
    }
    comp <- comp +
      geom_line(data = data.median, aes(x = Time, y = Frequency), color = color.plot, linewidth = 0.5) +
      geom_line(data = data.median, aes(x = Time, y = predicted_values.median), color = "darkred", linewidth = 0.2, alpha = 0.6) +
      geom_vline(xintercept = as.numeric(unlist(breakpoint.median)), linetype = 2, color = "black", linewidth = 0.2) +
      scale_x_continuous(expand = c(0,0), limits = c(0,20000), breaks = seq(0, 18000, 9000), name = "Years ago") +
      scale_y_continuous(expand = c(0,0), limits = c(0,1), name = "Frequency") +
      theme_minimal() +
      theme(aspect.ratio = 1,
            #axis.line = element_line(colour = "black", size = 0.2),
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            text = element_text(size = 7.5),
            legend.title = element_blank(),
            legend.text = element_text(size = 7.5),
            legend.spacing.y = unit(-0.2, 'cm'),
            axis.text = element_text(color = "black", size = 7.5),
            axis.title.x = element_text(color = "black", size = 7.5),
            axis.title.y = element_text(color = "black", size = 7.5),
            legend.key.size = unit(0.3, 'cm'),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            panel.grid = element_blank()
      )
    tiff_out = paste0(path, "/broken_stick_", unique(locus.name),"_npsi_", npsi.num,".test.tiff")
    tiff(tiff_out, units = "cm", res = 600, width = 8, height = 4)
    print(comp)
    dev.off()
    output <- c(locus = unique(locus.name), npsi = npsi.num, breakpoint = paste0(breakpoint.median, collapse = ','), breakpoint.min = paste0(breakpoint.min, collapse = ','), breakpoint.max = paste0(breakpoint.max, collapse = ','))
    output.txt <- rbind(output.txt, output)
    colnames(output.txt) <- c("locus", "npsi", "breakpoint", "breakpoint.min", "breakpoint.max")
  }
}
save(output.txt, file = paste0(path, "/broken_stick_results.rda"))
write.table(output.txt, file = paste0(path, "/broken_stick_results.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

