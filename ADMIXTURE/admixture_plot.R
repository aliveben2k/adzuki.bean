#args[1]: table file
#args[2]: output folder path
#args[3]: K value
library(reshape2)
library(ggplot2)
col.brewer.theme <- c("dodgerblue","orange","red","navy", "red4", "gold", "darkviolet", "maroon1", "peru", "coral1")
args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], sep="\t", header=TRUE)
K <- paste0('K=',as.character(args[3]))
col.pal <- col.brewer.theme[1:as.numeric(args[3])] #user_defined
df <- melt(df, id = "X")
a <- ggplot(df,aes(x=reorder(as.factor(X), as.numeric(variable)),y=value,fill=variable)) + 
  geom_bar(stat = "identity", width = 1) + 
  xlab(K) + 
  ylab(paste0("Value")) + 
  scale_x_discrete(labels=NULL, breaks=NULL) + 
  scale_fill_manual(values=col.pal) + 
  theme_light() + 
  theme(text = element_text(size=16), 
        panel.border = element_rect(colour = "white", fill=NA, size=0.5))
pdf(file = paste0(args[2],"/admixture_K",args[3],".pdf"), width = 12, height = 3)
print(a)
dev.off()