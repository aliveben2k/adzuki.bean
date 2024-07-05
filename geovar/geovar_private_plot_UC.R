library(ggplot2)
#Red series:
C_color = "navy"
R_color = "dodgerblue"
#Blue series: "darkblue","lightskyblue"
#Green series: "darkgreen","greenyellow"

args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], header = T, sep = "\t")
subset1 <- subset(data, grepl("[RC]", data$Group))
#subset1$X <- factor(subset1$X, levels = c("JP_LR", "CN_LRN", "CN_LRS", "JP_WL", "CN_WL", "SOUTH_WL2","SOUTH_WL1"))
subset1 <- subset1[order(subset1$X),]
subset1 <- t(subset1)
outtable <- sub("txt$", "C_only.txt", args[1])
write.table(subset1, file = outtable, sep = "\t", quote = F, col.names = F)
# geo.plot <- ggplot(subset1, aes(X, FV, fill = Group)) + 
#   geom_bar(stat="identity") + 
#   scale_fill_manual(values=c(C_color,R_color)) + 
#   theme_classic() + 
#   labs(x = 'Genogroup', y = 'Cumulative fraction of variants (%)', fill = 'Frequency') + 
#   theme(legend.position="right",
#         axis.title = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   scale_y_continuous(expand = c(0,0))
# output <- sub("txt$", "C_only.pdf", args[1])
# pdf(output, width=6, height=4)
# print(geo.plot)
# dev.off()