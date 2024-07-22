pckgs <- c("ggplot2","ggnewscale","RColorBrewer");
lapply(pckgs, library, character.only = TRUE)

treemix.bootstrap.m <- function(in.file, out.file = "tmp", phylip.file, pop.color.file = NULL, 
          nboot, label.size = 3, label.nudge = 0.001, x.expand = 0.01, flip = vector(), 
          arrow.size = 0.02, scale.ybar = 0.1, scale.xbar = 0, xmin = 0, line.size = 1, 
          scale = T, mbar = T, plotmig = T, plotnames = T, plotboot = T, text.size = 14,
          legend.location = NULL, legend.size = 12, mw.theme = NULL, node.size = 4, yaxis.show = F,
          ...)
{
  d = paste(in.file, ".vertices.gz", sep = "")
  e = paste(in.file, ".edges.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  e = read.table(gzfile(e), as.is = T, comment.char = "", quote = "")
  e[, 3] = e[, 3] * e[, 4]
  e[, 3] = e[, 3] * e[, 4]
  pops <- sort(d$V2[!is.na(d$V2)])
  if (length(pop.color.file) == 0) { #define population colors (o table)
    o <- as.data.frame(cbind(pops, rep("black", length(pops))), stringsAsFactors = F)
    colnames(o) <- c("V1", "V2")
  } else {
    if (!file.exists(pop.color.file)) {
      cat("File with population order and colors not found. Please check! \n")
      stop("Exit", call. = F)
    } else {
      o = read.table(pop.color.file, as.is = T, comment.char = "", quote = "")
      if (!setequal(x = pops, y = o[, 1])) {
        cat("File with population order and colors doesn't match with the population in the dataset. Please check! \n")
        stop("Exit", call. = F)
      }
    }
  }
  se = paste(in.file, ".covse.gz", sep = "")
  se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
  m1 = apply(se, 1, mean)
  mse = mean(m1)
  for (i in 1:length(flip)) {
    d = flip_node(d, flip[i])
  }
  d$x = "NA"
  d$y = "NA"
  d$ymin = "NA"
  d$ymax = "NA"
  options(warn = -1)
  d$x = as.numeric(d$x)
  d$y = as.numeric(d$y)
  d$ymin = as.numeric(d$ymin)
  d$ymax = as.numeric(d$ymax)
  options(warn = 0)
  d = set_y_coords(d)
  d = set_x_coords(d, e)
  tmplog <- capture.output({
    d = set_mig_coords(d, e)
  })
  vertices_treemix.df <- as.data.frame(matrix(NA, ncol = 1, nrow = nrow(d)))
  for (pos1 in 1:nrow(d)) {
    tmp1 <- unlist(strsplit(d[pos1, 11], split = ","))
    tmp3 <- NULL
    for (pos2 in 1:length(tmp1)) {
      tmp2 <- gsub(pattern = "\\(", replacement = "", x = unlist(strsplit(tmp1[pos2], ":"))[1])
      tmp3 <- c(tmp3, tmp2)
      rm(tmp2)
    }
    vertices_treemix.df[pos1, 1] <- paste(tmp3, collapse = " ")
    rm(tmp3)
  }
  vertices_treemix.df[vertices_treemix.df == "NA"] <- NA
  vertices_boot.df <- newick.split(phylip.file)
  d2 <- cbind(d, vertices_treemix.df)
  colnames(d2)[ncol(d2)] <- "treemix_vert"
  d2$Boot <- NA
  rstvert <- NULL
  for (b in 1:nrow(vertices_boot.df)) {
    vertsplit <- unlist(strsplit(vertices_boot.df[b, 2], split = " "))
    vertval <- vertices_boot.df[b, 1]
    vertlength <- length(vertsplit)
    pp <- 0
    for (a in vertsplit) {
      if (which(vertsplit == a) == 1) {
        tmpgrep <- d2$treemix_vert
      }
      tmpgrep2 <- tmpgrep[grep(pattern = a, x = tmpgrep)]
      tmpgrep <- tmpgrep2
      rm(tmpgrep2)
    }
    for (a in 1:length(tmpgrep)) {
      tmplng <- length(unlist(strsplit(tmpgrep[a], split = " ")))
      if (tmplng == vertlength) {
        d2$Boot[which(d2$treemix_vert == tmpgrep[a])] <- vertval
        pp <- 1
      }
    }
    if (pp == 0) {
      rstvert <- c(rstvert, vertices_boot.df[b, 2])
    }
    rm(tmpgrep)
  }
  if (length(rstvert) >= 1) {
    rstvert.df <- as.data.frame(rstvert)
    colnames(rstvert.df) <- "NoMatch"
    rstvert.df <- merge(x = rstvert.df, y = vertices_boot.df, by.x = "NoMatch", by.y = "POP", all.x = T, sort = F)
    rstvert.df$TreemixVertex <- NA
    notassigned <- d2$treemix_vert[which(is.na(d2$Boot))]
    if (length(notassigned) >= 1) {
      notassigned <- as.data.frame(notassigned)
      colnames(notassigned) <- "POP"
      notassigned$POP <- as.character(notassigned$POP)
      notassigned$LEN <- NA
      for (b in 1:nrow(notassigned)) {
        notassigned$LEN[b] <- length(unlist(strsplit(notassigned$POP[b], split = " ")))
      }
      notassigned <- notassigned[which(notassigned$LEN != 1), ]
    }
    for (b in rstvert.df$NoMatch) {
      vertsplit <- unlist(strsplit(b, split = " "))
      vertval <- vertices_boot.df[which(vertices_boot.df$POP == b), 1]
      vertlength <- length(vertsplit)
      check <- which(notassigned$LEN %in% c(vertlength))
      if (length(check) == 1) {
        treemix.vertsplit <- unlist(strsplit(notassigned$POP[check], split = " "))
        treemix.diff <- setdiff(treemix.vertsplit, vertsplit)
        tmp.vert <- c(vertsplit, treemix.diff)
        vertlength <- length(tmp.vert)
        for (a in tmp.vert) {
          if (which(tmp.vert == a) == 1) {
            tmpgrep <- vertices_boot.df$POP
          }
          tmpgrep2 <- tmpgrep[grep(pattern = a, x = tmpgrep)]
          tmpgrep <- tmpgrep2
          rm(tmpgrep2)
        }
        for (a in 1:length(tmpgrep)) {
          tmplng <- length(unlist(strsplit(tmpgrep[a], split = " ")))
          if (tmplng == vertlength) {
            rstvert.df$TreemixVertex[which(rstvert.df$NoMatch == b)] <- notassigned$POP[check]
            d2$Boot[which(d2$treemix_vert == notassigned$POP[check])] <- rstvert.df$VALUE[which(rstvert.df$NoMatch == b)]
          }
        }
      }
    }
    write.table(x = rstvert.df, file = paste(out.file, "_NoMatch_bootstrap.txt", sep = ""), quote = F, row.names = F, col.names = T)
  }
  if (!exists("nboot") & !is.numeric(nboot)) {
    cat("Number of bootstrap is unknown or is not a number. Please check! \n")
    stop("Exit", call. = F)
  }
  d2$Boot <- round(d2$Boot/nboot * 100, 2)
  d2$color <- NA
  d2$color[which(d2$Boot <= 100)] <- "a"
  d2$color[which(d2$Boot < 90)] <- "b"
  d2$color[which(d2$Boot < 75)] <- "c"
  d2$color[which(d2$Boot < 50)] <- "d"
  d2$point <- NA
  d2$point[which(d2$Boot <= 100)] <- 19
  d2$point[which(d2$Boot < 90)] <- 15
  d2$point[which(d2$Boot < 75)] <- 17
  d2$point[which(d2$Boot < 50)] <- 4
  d2$point <- as.numeric(d2$point)
  write.table(x = d2[which(d2[, 3] != "ROOT" & d2[, 4] != "MIG" & d2[, 5] != "TIP" & is.na(d2$Boot)), 16], 
              file = paste(out.file, "_NoMatch_treemix.txt", sep = ""), 
              quote = F, row.names = F, col.names = F)
  ##ggplot
  plot.ratio <- (max(d2$x)+x.expand-xmin)/max(d2$y)
  gg.plot <- ggplot() + labs(x = "Drift parameter", y = "") +
    scale_x_continuous(limits = c(xmin, max(d2$x)+x.expand)) +
    coord_fixed(ratio=plot.ratio)
  mw = max(e[e[, 5] == "MIG", 4]) #max weight
  if (mw > 0.5){
    scale.bar.max <- 1
    scale.bar.middle.name <- "0.5"
    scale.bar.max.name <- "1"
  } else {
    scale.bar.max <- 0.5
    scale.bar.middle.name <- "0.25"
    scale.bar.max.name <- "0.5"
  }
  mcols = rev(heat.colors(150)) #bar colors
  seg.coord.mit <- c()
  seg.coord.nomit <- c()
  for (i in 1:nrow(e)) { #plot tree line
    col = 1
    if (e[i, 5] == "MIG") {
      col = e[i, 4]
      if (is.na(col)) {
        col = NA
      }
    }
    v1 = d[d[, 1] == e[i, 1], ]
    v2 = d[d[, 1] == e[i, 2], ]
    if (e[i, 5] == "MIG") {
      if (plotmig) {
        seg.coord.mit.tmp <- cbind(v1[1, ]$x, v2[1,]$x, v1[1,]$y, v2[1,]$y, col)
        seg.coord.mit <- rbind(seg.coord.mit, seg.coord.mit.tmp)
      }
    } else {
      seg.coord.nomit.tmp <- cbind(v1[1, ]$x, v2[1,]$x, v1[1,]$y, v2[1,]$y, col)
      seg.coord.nomit <- rbind(seg.coord.nomit, seg.coord.nomit.tmp)
    }
  }
  colnames(seg.coord.mit) <- c("xstart","xend","ystart","yend","col")
  colnames(seg.coord.nomit) <- c("xstart","xend","ystart","yend","col")
  if (plotmig) {
    if (mbar){
      gg.plot <- gg.plot +
        geom_segment(data.frame(seg.coord.mit), arrow = arrow(length=unit(arrow.size,"npc"), type = "closed"), #admix group
                     mapping = aes(
                       x = as.numeric(xstart), y = as.numeric(ystart),
                       xend = as.numeric(xend), yend = as.numeric(yend),
                       color = col), 
                     size = line.size)
    } else {
      gg.plot <- gg.plot +
        geom_segment(data.frame(seg.coord.mit), arrow = arrow(length=unit(arrow.size,"npc"), type = "closed"), #admix group
                     mapping = aes(
                       x = as.numeric(xstart), y = as.numeric(ystart),
                       xend = as.numeric(xend), yend = as.numeric(yend),
                       color = col), 
                     size = line.size,
                     show.legend = F)
    }
  }
  if (grepl("rainbow", mw.theme)){
    if(!require("viridisLite")) install.packages("viridisLite")
    library(viridisLite)
    gg.plot <- gg.plot +
      scale_color_viridis_c(na.value="tan1",
                            option = "H",
                            direction = 1,
                            limits=c(0,scale.bar.max),
                            breaks=c(0,scale.bar.max/2,scale.bar.max),
                            labels=c("0", scale.bar.middle.name, scale.bar.max.name),
                            name = "Migration weight",
                            guide=guide_colourbar(ticks = F))
  } else {
    gg.plot <- gg.plot +
      scale_color_gradient(low="gold2",
                           high="purple1",
                           na.value="tan1",
                           limits=c(0,scale.bar.max),
                           breaks=c(0,scale.bar.max/2,scale.bar.max),
                           labels=c("0", scale.bar.middle.name, scale.bar.max.name),
                           name = "Migration weight",
                           guide=guide_colourbar(ticks = F))    
  }
  gg.plot <- gg.plot +
    geom_segment(data.frame(seg.coord.nomit), #pure group
                 mapping = aes(x = as.numeric(xstart), y = as.numeric(ystart),
                           xend = as.numeric(xend), yend = as.numeric(yend)),
                 size = line.size,
                 color = "black", 
                 show.legend = F)
  if (plotboot) { #plot bootstrep dots
    d2.point <- d2[!is.na(d2$color),]
    palette <- brewer.pal(9, "Blues")
    palette <- rev(palette[5:8])
    gg.plot <- gg.plot + new_scale_color() +
      geom_point(d2.point, mapping = aes(x = x, y = y, color = color, shape = factor(point)), size = node.size, na.rm = TRUE) +
      scale_shape_manual(values=c(19,15,17,4),
                         #na.value = NA, 
                         labels = c("90-100%","75-90%","50-75%","<50%"),
                         limits = factor(c(19,15,17,4)),
                         name = "Bootstrap values") +
      scale_color_manual(values=palette,
                         #na.value = NA, 
                         labels = c("90-100%","75-90%","50-75%","<50%"),
                         limits = c("a","b","c","d"),
                         name = "Bootstrap values")
  }
  tmp = d[d[, 5] == "TIP", ]
  pop.labels <- c()
  for (i in 1:nrow(tmp)) { #put population lebels on the nodes
    tcol = o[o[, 1] == tmp[i, 2], 2]
    if (plotnames) {
      pop.labels.tmp <- cbind(tmp[i, ]$x, tmp[i, ]$y, tmp[i,2], tcol)
      pop.labels <- rbind(pop.labels, pop.labels.tmp)
    }
  }
  if (length(pop.labels) != 0){
    colnames(pop.labels) <- c("x.pos","y.pos","labs","col")
    gg.plot <- gg.plot + new_scale_color() +
      geom_text(data.frame(pop.labels),
                mapping = aes(x = as.numeric(x.pos), 
                              y = as.numeric(y.pos),
                              color = col,
                              label = labs),
                show.legend = F,
                nudge_x = label.nudge,
                hjust = 0,
                size = label.size) +
      scale_color_manual(values=pop.labels[,4],
                         breaks=pop.labels[,4],
                         na.value = NA)
  }
  if (scale) { #label.nudgelay s.e. bar
    gg.plot <- gg.plot +
      geom_segment(mapping = aes(x = scale.xbar, 
                                 xend = scale.xbar + (mse * 10), 
                                 y = scale.ybar, yend = scale.ybar),
                   size = line.size) +
      geom_text(mapping = aes(x = scale.xbar + (mse * 10)/2, 
                              y = scale.ybar, label = "10 s.e."), 
                show.legend = F, 
                nudge_y = -0.03, 
                size = label.size, 
                hjust = 0.5)
  }
  if (legend.size < 4){
    legend.size <- 4
  }
  gg.plot <- gg.plot +
    theme_classic() +
    theme(text = element_text(size = text.size),
          axis.text.x = element_text(color = "black"),
          legend.title = element_text(face = "bold", size = legend.size),
          legend.text = element_text(size = legend.size - 1))
  if (yaxis.show == F){
    gg.plot <- gg.plot +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank())
  }
  if (length(legend.location) != 0){
    if (legend.location == "bottomright"){
      gg.plot <- gg.plot +
        theme(
          legend.justification = c(1, 0), 
          legend.position = c(1, 0)
        )
    } else if (legend.location == "topright"){
      gg.plot <- gg.plot +
        theme(
          legend.justification = c(1, 1), 
          legend.position = c(1, 1)
        )
    } else if (legend.location == "bottomleft"){
      gg.plot <- gg.plot +
        theme(
          legend.justification = c(0, 0), 
          legend.position = c(0, 0)
        )
    } else if (legend.location == "topleft"){
      gg.plot <- gg.plot +
        theme(
          legend.justification = c(0, 1), 
          legend.position = c(0, 1)
        )
    } else {
      gg.plot <- gg.plot +
        theme(legend.position = legend.location)
    }
  }
  print(gg.plot)
}
