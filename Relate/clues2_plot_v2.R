library(argparse)
library(ggplot2)
library(scales)
library(gridExtra)
library(vroom)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument('--freqs', type="character")
parser$add_argument('--post', type="character")
parser$add_argument('--figure', type="character")
parser$add_argument('--posterior_intervals', type="double", nargs='+', 
                    default=c(0.5, 0.75, 0.95, 0.999),
                    help='The posterior thresholds for which to draw different consensus trees.')
parser$add_argument('--generation_time', type="double", default=-1.0, 
                    help='Conversion to generation time.')
args <- parser$parse_args()

# Prepare confidence intervals and color intensities
ConfidenceIntervals <- sort(round(args$posterior_intervals, 4))
ColorIntensity <- seq_along(ConfidenceIntervals) / length(ConfidenceIntervals)

# Load data
freqs <- as.numeric(unlist(read.csv(args$freqs, header = FALSE)))
logpost <- as.matrix(vroom(args$post, delim = ",", col_names = FALSE, col_types = "d")) #use vroom to speed up the reading
epochs <- seq(0, ncol(logpost))
colnames(logpost) <- epochs[1:length(epochs)-1]
y.max <- apply(logpost, 2, which.max) #find the index of the maximum value in each epochs
y.max <- freqs[y.max] #get real frequency of the epochs

# Initialize variables
MATRIXTOPLOT <- matrix(0, nrow = nrow(logpost), ncol = ncol(logpost))
StartIndices <- rep(1, length(epochs))

# Determine start indices
for (timeinterval in head(epochs, -1)) {
  for (lowerfrequencyindex in seq_len(nrow(logpost))) {
    if (logpost[lowerfrequencyindex, timeinterval + 1] > 0.000001 / length(freqs)) {
      StartIndices[timeinterval + 1] <- max(1, lowerfrequencyindex)
      break
    }
  }
}

# Calculate posterior intervals
for (indexx in seq_along(ConfidenceIntervals)) {
  for (timeinterval in head(epochs, -1)) {
    tentativespan <- Inf
    currentlow <- -1
    currenthigh <- -1
    TimeSLICE <- logpost[, timeinterval + 1]
    
    for (lowerfrequencyindex in seq(StartIndices[timeinterval + 1], nrow(logpost)-1)) {
      possiblesum <- TimeSLICE[lowerfrequencyindex]
      
      for (higherfrequencyindex in (lowerfrequencyindex + 1):nrow(logpost)) {
        possiblesum <- possiblesum + TimeSLICE[higherfrequencyindex]
        rangeofsum <- freqs[higherfrequencyindex] - freqs[lowerfrequencyindex]
        
        if (possiblesum >= ConfidenceIntervals[indexx] && rangeofsum < tentativespan) {
          tentativespan <- rangeofsum
          currentlow <- lowerfrequencyindex
          currenthigh <- higherfrequencyindex
        } else if (rangeofsum > tentativespan) {
          break
        }
      }
    }
    
    for (j in seq_len(nrow(MATRIXTOPLOT))) {
      if (j >= currentlow && j <= currenthigh) {
        MATRIXTOPLOT[j, timeinterval + 1] <- max(MATRIXTOPLOT[j, timeinterval + 1], ColorIntensity[indexx])
      }
    }
  }
}

# # Convert MATRIXTOPLOT to a data frame
plot_data <- melt(MATRIXTOPLOT)
colnames(plot_data) <- c("freq_index", "epoch_index", "intensity")
plot_data$freq <- freqs[plot_data$freq_index] #add real frequency value according to the frequency index
plot_data$epoch <- (plot_data$epoch_index - 1) * ifelse(args$generation_time == -1.0, 1, args$generation_time)

# Plot using ggplot2
ggplot(plot_data, aes(x = epoch, y = freq, fill = intensity)) +
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(100, "YlOrRd"), name = "Intensity") +
  labs(
    x = ifelse(args$generation_time == -1.0, "Generations before present", "Years before present"),
    y = "Derived allele frequency",
    title = "Posterior Density Heatmap"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )









rownames(MATRIXTOPLOT) <- freqs
plot.x <- c(); plot.y <- c()
for (i in 1:ncol(MATRIXTOPLOT)){
  maxies <- c()
  for (j in 1:nrow(MATRIXTOPLOT)){
   if (MATRIXTOPLOT[j,i] == 1){
     maxies <- c(maxies, j)
   }
  }
  maxies <- sort(maxies, decreasing = F)
  maxies.freqs <- as.numeric(rownames(MATRIXTOPLOT)[maxies])
  median.freq <- median(maxies.freqs, na.rm = T)
  plot.x <- c(plot.x, epochs[i])
  plot.y <- c(plot.y, median.freq)
}
if (args$generation_time != -1.0) {
  plot.x <- plot.x * args$generation_time
  xlab = "Years before present"
} else {
  xlab = "Generations before present"
}
y.max <- cbind(x = plot.x, y = plot.y)

# Prepare data for plotting
# plot_data <- melt(MATRIXTOPLOT)
# freqs <- rep(freqs, length(epochs) - 1)
# plot_data$Var1 <- freqs
# colnames(plot_data) <- c("Frequency", "Time", "value")
#plot_data <- expand.grid(freqs = freqs, epochs = epochs[-length(epochs)])
#plot_data$intensity <- as.vector(MATRIXTOPLOT)

# Plot results
# if (args$generation_time == -1.0) {
#   image(1:(length(epochs) - 1), freqs, t(MATRIXTOPLOT), col = hcl.colors(100, "YlOrRd"), xlab = "Generations before present",
#         ylab = "Derived allele frequency", useRaster = TRUE)
# } else {
#   image(1:(length(epochs) - 1) * args$generation_time, freqs, t(MATRIXTOPLOT), col = hcl.colors(100, "YlOrRd"), 
#         xlab = "Years before present", ylab = "Derived allele frequency")
# }

# # Plot
# if (args$generation_time == -1.0) {
#   x_label <- "Generations before present"
#   plot_data$Time <- plot_data$Time
# } else {
#   x_label <- "Years before present"
#   plot_data$Time <- plot_data$Time * args$generation_time
# }

# p <- ggplot(plot_data, aes(x = Time, y = freqs, fill = intensity)) +
#   geom_tile() +
#   scale_fill_gradientn(colors = rev(rainbow(length(ColorIntensity), end = 0.7)), 
#                        #values = rescale(ColorIntensity),
#                        name = "Posterior Interval") +
#   theme_minimal() +
#   labs(x = x_label, y = "Derived allele frequency") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12))

p <- ggplot(post.exp, aes(x=Time, y=Frequency, fill=value)) + 
  geom_tile() +
  #geom_line(data = y.max, aes(x = x, y = as.numeric(y)), color = "black", inherit.aes = FALSE, linewidth = 0.2) +
  #geom_line(data = y.max, aes(x = x, y = as.numeric(y)), color = "black", inherit.aes = FALSE, size = 0.2) +
  #geom_smooth(data = y.max, aes(x = x, y = as.numeric(y)), color = "black", method = "gam", inherit.aes = FALSE, linewidth = 0.5) +
  scale_fill_distiller(palette = color.series, 
                       direction = 1,
                       name = "Posterior probability",
                       guide = guide_colourbar(title.position = "right", 
                                               barheight = unit(fig_height*0.8, "inches"), 
                                               barwidth = unit(fig_width*0.01, "inches"),
                                               ticks = F,
                                               ticks.color = "black",
                                               ticks.linewidth = 1,
                                               label.vjust = 0.5, 
                                               title.hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        text = element_text(face = "bold", size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0)),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0), angle = -90),
        legend.text = element_text(size = round(4*min(c(fig_width, fig_height), na.rm = T), digits = 0)),
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5))



# Save plot
ggsave(filename = paste0(args$figure, ".png"), plot = p, dpi = 300, width = 20, height = 10)