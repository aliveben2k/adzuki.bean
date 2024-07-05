#only for Japan map
options(warn=-1)
s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Jul. 2022
Usage: Rscript draw_genetic_median_dist.R -l INFO_LIST_FILE -m DISTANCE_MATRIX -t POP_NAME [-s POP_NAME] [-md RESULT_FILE] [-r COORDINATE] [-fs WIDTH_HEIGHT] [-tm THEME] [-lp POSITION] [-ts TEXT_SIZE] [-k NUMBER] [-ktm THEME] [-cty NAME] [-kex EXCLUDE_SUBREGION] [--kmodel MODEL] [--kpsill NUMBER] [--knugget NUMBER] [--krange NUMBER] [--kregion MAP_BORDER] [-h]\n
-l/--list: Information list file. Table seperates by tab with header line. format: ID Population Long. Lat.
-m/--matrix: genetic distance matrix file.
-t/--target: target population.
-s/--source: source population.
-ts/--text_size: text size.
-md/--mdist_file: result table of this analysis. If provided, -l and -m are not necessary.
-r/--region: map coordinate (left,right,bottom,top). default: -180,180,-90,90
-fs/--figure_size: figure size (width,height). default: 10,6
-tm/--theme: map theme. Possible values: light, dark. default: light
-lp/--legend_position: legend position. Possible values: topright, bottomright, topleft, bottomleft. default: system default
-k/--kriging_res: Kriging resolution in kilometer. default: 10. Be careful: smaller number may take more time.
-cty/--country: country name that you are looking at. Multiple regions can be seperated by comma.
-kex/--kexclude_subregion: exclude subregions from analysis (Japan only), kriging only. Multiple subregions can be seperated by comma.
-ktm/--ktheme: map theme. Possible values: rainbow, light, dark. default: rainbow
--kmodel: model used for Krigging. default: Sph
--kpsill: psill value used for Krigging. default: interval formula
--knugget: nugget value used for Krigging. default: interval formula
--krange: range value used for Krigging. default: interval formula
--kregion: map coordinate border (left,right,bottom,top) for kriging plot. default: -180,180,-90,90
  **If you got incorrect output, increasing border range of --kregion might solve the problem.
-h/--help: help.\n")
}
get.target.idx <- function(x, y){
  f.target.idx <- c()
  for (m in 1:length(y)){
    if (y[m] == x){
      f.target.idx <- m
    }
  }
  if (length(f.target.idx) == 0){
    cat("Cannot find the target population in the list.\n")
    quit()
  }
  return(f.target.idx)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
mdist.table <- c()
map.theme <- "light"
text.size <- c()
kg <- c()
krig.map.theme <- "rainbow"
cty <- c()
exclude.region <- c()
krig.range.in <- c()
krig.nugget.in <- c()
krig.psill.in <- c()
krig.model <- c()
krig.region <- c()
region <- c()
fig_size <- c()
target <- c()
source <- c()
leg.pos <- c()
leg.pos.values <- c("bottomright", "bottomleft", "topright", "topleft")
for (i in 1:length(args)){
  if (args[i] == '-l' || args[i] == '--list'){
    info <- args[i+1]
    if (!file.exists(info)){
      s.help()
      cat("-l: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-m' || args[i] == '--matrix'){
    dist.matrix <- args[i+1]
    if (!file.exists(dist.matrix)){
      s.help()
      cat("-m: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-h' || args[i] == '--help'){
    s.help()
    quit()
  }
  if (args[i] == '-t' || args[i] == '--target'){
    target <- args[i+1]
  }
  if (args[i] == '-s' || args[i] == '--source'){
    source <- args[i+1]
  }
  if (args[i] == '-ts' || args[i] == '--text_size'){
    text.size <- as.numeric(args[i+1])
  }
  if (args[i] == '-kex' || args[i] == '--kexclude_subregion'){
    exclude.region <- args[i+1]
  }
  if (args[i] == '-r' || args[i] == '--region'){
    region <- args[i+1]
    #'128,150,30,46' for japan
  }
  if (args[i] == '--kregion'){
    krig.region <- args[i+1]
    #128,150,30,46 for japan
    #70,150,15,57 for V. angularis
  }
  if (args[i] == '-fs' || args[i] == '--figure_size'){
    fig_size <- args[i+1]
    #'6.5,4.7' for japan
  }
  if (args[i] == '-k' || args[i] == '--kriging_res'){
    if (!grepl("^-", args[i+1]) && i < length(args)){
      kg <- as.numeric(args[i+1])
      #you can try 10 first
      kg <- kg * 1000
    } else {
      kg <- 10000
    }
  }
  if (args[i] == '-cty' || args[i] == '--country'){
    cty <- as.character(args[i+1])
  }
  if (args[i] == '--kpsill'){
    krig.psill.in <- as.numeric(args[i+1])
  }
  if (args[i] == '--knugget'){
    krig.nugget.in <- as.numeric(args[i+1])
  }
  if (args[i] == '--krange'){
    krig.range.in <- as.numeric(args[i+1])
  }
  if (args[i] == '--kmodel'){
    krig.model <- as.character(args[i+1])
  }
  if (args[i] == '-md' || args[i] == '--m_dist_file'){
    mdist.table <- as.character(args[i+1])
    if (!file.exists(mdist.table)){
      s.help()
      cat("-md: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-tm' || args[i] == '--theme'){
    if (grepl("light|dark", args[i+1])){
      map.theme <- args[i+1]
    } else {
      cat("-tm: cannot recognize the theme used. Default setting used.\n")
    }
  }
  if (args[i] == '-lp' || args[i] == '--legend_position'){
    leg.pos.in <- args[i+1]
    for (j in 1:length(leg.pos.values)){
      if (leg.pos.in == leg.pos.values[j]){
        leg.pos = leg.pos.in
      }
    }
    if (length(leg.pos) < 1){
      cat("-lp: cannot recognize the value. Default setting used.\n")
    }
  }
  if (args[i] == '-ktm' || args[i] == '--ktheme'){
    if (grepl("light|dark|rainbow", args[i+1])){
      krig.map.theme <- args[i+1]
    } else {
      cat("-ktm: cannot recognize the theme used. Default setting used.\n")
    }
  }
}
if (length(mdist.table) < 1){
  if (length(info) < 1){
    s.help()
    cat("-l: argument must be provided.\n")
    quit()
  }
  if (length(dist.matrix) < 1){
    s.help()
    cat("-m: argument must be provided.\n")
    quit()
  }
}
if (length(target) < 1){
  s.help()
  cat("-t: argument must be provided.\n")
  quit()
}
if (length(kg) > 0){
  if(!require("gstat")) install.packages("gstat")
  if(!require("sf")) install.packages("sf")
  library(gstat)
  library(sf)
  if (length(krig.model) == 0){
    krig.model <- "Sph"
  }
  #sf::sf_use_s2(FALSE)
}
if (length(krig.range.in) != 0 || length(krig.nugget.in) != 0 || length(krig.psill.in) != 0){
  if (length(kg) == 0){
    cat("-k: Kriging option is not set. --krange, --kpsill or --knugget will not be used.\n")
  }
}

#set region size
left <- -180
right <- 180
bottom <- -90
top <- 90
if (length(region) != 0){
  coord <- unlist(strsplit(region, ','))
  left <- as.numeric(coord[1])
  right <- as.numeric(coord[2])
  bottom <- as.numeric(coord[3])
  top <- as.numeric(coord[4])
}
krig.left <- c()
krig.right <- c()
krig.bottom <- c()
krig.top <- c()
if (length(krig.region) != 0){
  krig.coord <- unlist(strsplit(krig.region, ','))
  krig.left <- as.numeric(krig.coord[1])
  krig.right <- as.numeric(krig.coord[2])
  krig.bottom <- as.numeric(krig.coord[3])
  krig.top <- as.numeric(krig.coord[4])
}
#set figure size
fig_width <- 10
fig_height <- 7
if (length(fig_size) != 0){
  size <- unlist(strsplit(fig_size, ','))
  fig_width <- as.numeric(size[1])
  fig_height <- as.numeric(size[2])
}

#set text size
if (length(text.size) == 0){
  text.size <- round(3.5*min(c(fig_width, fig_height), na.rm = T), digits = 0)
}

if (length(cty) != 0){
  cty <- unlist(strsplit(cty, ','))
} else {
  cty <- '.'
}

if (length(source) != 0){
  source <- unlist(strsplit(source, ','))
}

#BEGIN of the script
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("geosphere")) install.packages("geosphere")
if(!require("sp")) install.packages("sp")
if(!require("sf")) install.packages("sf")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("viridisLite")) install.packages("viridisLite")
if(!require("mapdata")) install.packages("mapdata")
library(tidyverse)
library(geosphere)
library(sp)
library(sf)
library(ggplot2)
library(viridisLite)
library(mapdata)

if (length(mdist.table) < 1){
  data <- read.table(info, header = T, sep = "\t")
  pops <- unique(data[,2])
  long <- colnames(data)[3]
  lat <- colnames(data)[4]
  data.coord <- as.matrix(data[,3:4])
  colnames(data.coord) <- c(long, lat)
  rownames(data.coord) <- data[,1]
  geo.dist.matrix <- distm(data.coord, fun = distGeo) #get the geo distance matrix
  colnames(geo.dist.matrix) <- data[,1]
  rownames(geo.dist.matrix) <- data[,1]
  geo.dist.matrix <- as.data.frame(geo.dist.matrix)
  legend.cnt <- length(unique(data[,2]))
  
  gen.dist.matrix <- read.table(dist.matrix, header = T, sep = "\t", row.names = 1) #genetic distance matrix
  
  gen.matrix.index <- c()
  for (i in 1:nrow(data)){ #subset the matrix
    for (j in 1:nrow(gen.dist.matrix)){
      if (data[i,1] == rownames(gen.dist.matrix)[j]){
        gen.matrix.index <- c(gen.matrix.index, j)
      }
    }
  }
  #gen.matrix.index <- sort(gen.matrix.index)
  gen.dist.matrix <- gen.dist.matrix[gen.matrix.index, gen.matrix.index]
  
  #make a blank data frame
  out.table <- matrix(ncol = 6, nrow = 0)
  out.table <- matrix(ncol = 6, nrow = 0)
  out.table <- data.frame(out.table)
  
  #get the target population sample index
  target.idx <- get.target.idx(target, pops)
  target.sample.idx <- c()
  for (m in 1:nrow(data)){
    if (data[m,2] == target){
      target.sample.idx <- c(target.sample.idx, m)
    }
  }
  #get the target population sample index
  source.idx <- c()
  if (length(source) != 0){
    for (k in 1:length(source)){
      source.idx.tmp <- get.target.idx(source[k], pops)
      source.idx <- c(source.idx, source.idx.tmp)
    }
  }
  for (n in 1:length(pops)){
    if (n == target.idx){
      next
    }
    if (length(source) != 0){
      exist <- 0
      for (l in 1:length(source.idx)){
        if (n == source.idx[l]){
          exist <- 1
        }
      }
      if (exist == 0){
        next
      }
    }
    current.sample.idx <- c()
    for (p in 1:nrow(data)){ #get the current population sample index
      if (data[p,2] == pops[n]){
        current.sample.idx <- c(current.sample.idx, p)
      }
    }
    subset.gen.dist.matrix <- gen.dist.matrix[target.sample.idx, current.sample.idx] #get target-current genetice distance matrix
    for (k in 1:nrow(subset.gen.dist.matrix)){
      sample.median <- median(as.numeric(subset.gen.dist.matrix[k,]), na.rm = T)
      for (r in 1:nrow(data)){
        if (rownames(subset.gen.dist.matrix)[k] == data[r,1]){
          sample.coord <- c(data[r,3], data[r,4])
        }
      }
      temp <- c(rownames(subset.gen.dist.matrix)[k], target, pops[n], sample.median, sample.coord)
      out.table <- rbind(out.table, temp)
    }
  }
  colnames(out.table) <- c("ID", "sample_pop", "target_pop", "gen_m_dist", "long", "lat")
  out <- gsub("txt$", paste0(target, "_m_dist.table.txt"), info)
  write.table(out.table, file = out, quote = F, col.names = T, row.names = F, sep = "\t")
} else {
  out.table <- read.table(mdist.table, header = T, sep = "\t")
  if (length(colnames(out.table)) != 6){
    s.help()
    cat("-md: File format cannot be recognized.\n")
    quit()
  }
  pops <- unique(unlist(out.table[,2:3]))
  if (length(source) != 0){
    out.table.sub <- c()
    for (k in 1:length(source)){
      out.table.sub <- data.frame(rbind(out.table.sub, out.table[which(out.table$target_pop == source[k]),]))
    }
    out.table <- out.table.sub
    pops <- unique(unlist(out.table[,2:3]))
  }
  target.idx <- get.target.idx(target, pops)
  info <- mdist.table
  info <- gsub(paste0(target, "_m_dist.table.txt"), "txt", info)
}

for (n in 1:length(pops)){
  if (n == target.idx){
    next
  }
  world.sf <- st_as_sf(maps::map("world2", plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) #load world map
  subset.out.table <- subset(out.table, out.table[,3]==pops[n])
  subset.out.table[,4] <- as.numeric(subset.out.table[,4])
  subset.out.table[,5] <- as.numeric(subset.out.table[,5])
  subset.out.table[,6] <- as.numeric(subset.out.table[,6])
  long <- subset.out.table[,5]
  lat <- subset.out.table[,6]
  if (grepl("light", map.theme)){
    #light theme
    #out <- gsub("txt$", paste0(trait.name, ".light_theme.tiff"), file)
    out <- gsub("txt$", paste0(target, "_vs_",pops[n],".light_theme.tiff"), info)
    out <- gsub("\\.+",".", out)
    theme.country.line.color <- "black"
    theme.country.fill <- "white"
    theme.background <- "lightgray"
    theme.color.low <- "gold2"
    theme.color.high <- "blue"
  } else {
    #dark theme
    #out <- gsub("txt$", paste0(trait.name, ".dark_theme.tiff"), file)
    out <- gsub("txt$", paste0(target, "_vs_",pops[n],".dark_theme.tiff"), info)
    out <- gsub("\\.+",".", out)
    theme.country.line.color <- "grey50"
    theme.country.fill <- "black"
    theme.background <- "grey85"
    theme.color.low <- "gold2"
    theme.color.high <- "purple1"
  }
  out.plot <- ggplot() + 
    geom_sf(data = world.sf,
            color = theme.country.line.color, 
            fill = theme.country.fill, 
            size = 0.1) +
    xlim(left,right) + 
    ylim(bottom,top) +
    labs(title = paste0(target, " vs. ",pops[n])) +
    geom_point(data = subset.out.table, 
               aes(long, 
                   lat, 
                   color = subset.out.table[,4]), 
               alpha =0.9, 
               shape = 16, 
               size = round(0.1*min(c(fig_width, fig_height), na.rm = T), digits = 1)) +
    theme(text = element_text(size = text.size),
          panel.background = element_rect(fill = theme.background, color = theme.background),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
          axis.text = element_text(color = "black", size = text.size-2),
          axis.ticks = element_line(colour = "black", size = 0.2),
          axis.title = element_text(color = "black", size = text.size),
          #legend.title = element_text(size = text.size-2),
          #legend.title = element_blank(),
          legend.text = element_text(size = text.size),
          legend.key = element_rect(fill = NA, color = NA),
          legend.key.height= unit(fig_height*0.1, 'cm'),
          legend.key.width= unit(fig_width*0.05, 'cm'),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Longitude") + 
    ylab("Latitude")
  out.plot <- out.plot +
    scale_color_gradient(low=theme.color.low, 
                         high=theme.color.high, 
                         name = "Genetic\ndistance", 
                         #breaks = breaks, 
                         #labels = format(breaks)
    )
  if (length(leg.pos) > 0){
    if (leg.pos == "bottomright"){
      out.plot <- out.plot +
        theme(
          legend.justification = c(1, 0), 
          legend.position = c(1, 0)
        )
    }
    if (leg.pos == "topright"){
      out.plot <- out.plot +
        theme(
          legend.justification = c(1, 1), 
          legend.position = c(1, 1)
        )
    }
    if (leg.pos == "bottomleft"){
      out.plot <- out.plot +
        theme(
          legend.justification = c(0, 0), 
          legend.position = c(0, 0)
        )
    }
    if (leg.pos == "topleft"){
      out.plot <- out.plot +
        theme(
          legend.justification = c(0, 1), 
          legend.position = c(0, 1)
        )
    }
  }
  #pdf(out, width=fig_width, height=fig_height)
  tiff(out, units="cm", width = fig_width, height = fig_height, res = 600)
  print(out.plot)
  dev.off()
  
  if (length(kg) > 0){
    if (length(krig.left) == 0){
      krig.left <- left
      krig.right <- right
      krig.bottom <- bottom
      krig.top <- top
    }
    if (length(cty) != 0){
      cty <- unlist(strsplit(cty, ','))
    } else {
      cty <- '.'
    }
    subset.out.table.coord <- subset.out.table
    subset.out.table.coord$gen_m_dist <- subset.out.table.coord$gen_m_dist*100
    coordinates(subset.out.table.coord) <- ~long+lat
    #remove any duplicate coordinate
    subset.out.table.coord <- subset.out.table.coord[which(!duplicated(subset.out.table.coord@coords)),]
    #transform coordinate system to meter unit
    proj4string(subset.out.table.coord) <- CRS("+init=epsg:4326")
    CRS.new <- CRS("+init=epsg:3857")
    subset.out.table.coord.meter <- spTransform(subset.out.table.coord, CRS.new)
    #transform data to vgm variables
    subset.out.table.vgm <- variogram(gen_m_dist~1, subset.out.table.coord.meter)
    #subset.out.table.vgm <- variogram(subset.out.table.coord.meter@data[[gen_m_dist]]~1, subset.out.table.coord.meter)
    #remove any duplicate coordinate
    #subset.out.table.coord.meter <- subset.out.table.coord.meter[which(!duplicated(subset.out.table.coord.meter@coords)),]
    #set up Kriging parameters
    if (length(krig.range.in) == 0){
      krig.range <- round(max(subset.out.table.vgm$dist)-max(subset.out.table.vgm$dist)/100)
    } else {
      krig.range <- krig.range.in
    }
    krig.sill.digit <- ceiling(-log10(max(subset.out.table.vgm$gamma)))
    krig.sill <- quantile(subset.out.table.vgm$gamma, 0.9)
    if (length(krig.nugget.in) == 0){
      krig.nugget <- as.numeric(quantile(subset.out.table.vgm$gamma[1:5], 0.1))-((krig.sill-as.numeric(quantile(subset.out.table.vgm$gamma[1:5], 0.1)))/length(subset.out.table.vgm$gamma))
    } else {
      krig.nugget <- krig.nugget.in
    }
    if (length(krig.psill.in) == 0){
      if (krig.sill - krig.nugget > 0){
        krig.psill <- krig.sill - krig.nugget
      } else {
        cat("--knugget value is too big. Reseted by the interval formula.\n")
        krig.nugget <- as.numeric(quantile(subset.out.table.vgm$gamma[1:5], 0.1))-((krig.sill-as.numeric(quantile(subset.out.table.vgm$gamma[1:5], 0.1)))/length(subset.out.table.vgm$gamma))
        krig.psill <- krig.sill - krig.nugget
      }
    } else {
      krig.psill <- krig.psill.in
    }
    #apply model (Exp, Sph, Lin or Gau)
    model.krig <- vgm(krig.psill, model=krig.model, nugget=krig.nugget, range=krig.range)
    #fit model
    fit.variog <- fit.variogram(subset.out.table.vgm, model.krig, fit.kappa = T, warn.if.neg = F)
    #plot fitness of the model
    fitness <- plot(subset.out.table.vgm, fit.variog, name = "Genetic\ndistance")
    out <- gsub("txt$", paste0(target, "_vs_",pops[n],".k_model.fitness.tiff"), info)
    #pdf(out, width=4, height=3)
    tiff(out, units="cm", width = 12, height = 9, res = 600)
    print(fitness)
    dev.off()
    #set up extend region
    range.coord <- as.data.frame(c(krig.left,krig.right))
    range.coord <- cbind(range.coord, c(krig.bottom,krig.top))
    colnames(range.coord) <- c("long","lat")
    coordinates(range.coord) <- ~long+lat
    proj4string(range.coord) <- CRS("+init=epsg:4326")
    range.coord.meter <- spTransform(range.coord, CRS.new)
    range.coord.meter <- as.data.frame(range.coord.meter)
    x.range.min <- c()
    x.range.max <- c()
    y.range.min <- c()
    y.range.max <- c()
    if (min(range.coord.meter[,1]) - min(range.coord.meter[,1])*0.2 < -20037508){
      x.range.min <- -20037508
    } else {
      x.range.min <- min(range.coord.meter[,1]) - min(range.coord.meter[,1])*0.2
    }
    if (max(range.coord.meter[,1]) + max(range.coord.meter[,1])*0.2 > 20037508){
      x.range.max <- 20037508
    } else {
      x.range.max <- max(range.coord.meter[,1]) + max(range.coord.meter[,1])*0.2
    }
    if (min(range.coord.meter[,2]) - min(range.coord.meter[,2])*0.2 < -88985946){
      y.range.min <- -88985946
    } else {
      y.range.min <- min(range.coord.meter[,2]) - min(range.coord.meter[,2])*0.2
    }
    if (max(range.coord.meter[,2]) + max(range.coord.meter[,2])*0.2 > 88985946){
      y.range.max <- 88985946
    } else {
      y.range.max <- max(range.coord.meter[,2]) + max(range.coord.meter[,2])*0.2
    }
    x.range <- as.integer(c(min(range.coord.meter[,1]),max(range.coord.meter[,1])))
    y.range <- as.integer(c(min(range.coord.meter[,2]),max(range.coord.meter[,2])))
    grd <- expand.grid(long=seq(from=x.range[1], to=x.range[2], by=kg), lat=seq(from=y.range[1], to=y.range[2], by=kg))
    coordinates(grd) <- ~long+lat
    proj4string(grd) <- CRS("+init=epsg:3857")
    gridded(grd) <- TRUE
    #print(length(grd$lat))
    #do Kriging modelling
    cat("Applying Kriging model to the map...\nIt may take some time...\n")
    subset.out.table.kriged <- krige(gen_m_dist~1, subset.out.table.coord.meter, grd, model = model.krig)
    #subset.out.table.kriged <- krige(subset.out.table.coord.meter@data[[gen_m_dist]]~1, subset.out.table.coord.meter, grd, model = model.krig)
    proj4string(subset.out.table.kriged) <- CRS("+init=epsg:3857")
    subset.out.table.kriged <- spTransform(subset.out.table.kriged, CRS("+init=epsg:4326"))
    subset.out.table.kriged@data[["var1.pred"]] <- subset.out.table.kriged@data[["var1.pred"]]/100
    subset.out.table.kriged@data[["var1.var"]] <- subset.out.table.kriged@data[["var1.var"]]/100
    subset.out.table.kriged <- as.data.frame(subset.out.table.kriged)
    #do masking of the rest region except target region
    sf_use_s2(FALSE)
    region.map.mask <- st_as_sf(maps::map("japan", plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) %>% st_geometry() %>% st_cast('POLYGON')
    region.map.mask <- st_make_valid(region.map.mask)
    region.map.mask <- region.map.mask %>% st_union()
    #region.map.mask <- region.map.mask %>% st_set_crs(4326) %>% st_geometry() %>% st_cast('POLYGON') %>% st_union()
    region.bbox <- as.data.frame(t(bbox(grd)))
    region.bbox[1,1] <- region.bbox[1,1] - region.bbox[1,1]*0.2
    if (region.bbox[1,1] < -20037508){
      region.bbox[1,1] <- -20037508
    }
    region.bbox[2,1] <- region.bbox[2,1] + region.bbox[2,1]*0.2
    if (region.bbox[2,1] > 20037508){
      region.bbox[2,1] <- 20037508
    }
    region.bbox[1,2] <- region.bbox[1,2] - region.bbox[1,2]*0.2
    if (region.bbox[1,2] < -88985946){
      region.bbox[1,2] <- -88985946
    }
    region.bbox[2,2] <- region.bbox[2,2] + region.bbox[2,2]*0.2
    if (region.bbox[2,2] > 88985946){
      region.bbox[2,2] <- 88985946
    }
    region.bbox <- st_as_sf(region.bbox, coords = c(1:2)) %>% st_set_crs(3857) #meter coord system
    region.bbox <- st_transform(region.bbox, crs = 4326)
    region.bbox <- st_bbox(region.bbox) %>% st_as_sfc() %>% st_as_sf()
    region.map.mask <- st_difference(region.bbox, region.map.mask)
    
    #world.sf.line <- world.sf
    for (q in 1:length(cty)){
      world.sf.num.before <- nrow(world.sf)
      world.sf <- world.sf[!grepl(cty[q], world.sf$ID, ignore.case = T),]
      world.sf.num.after <- nrow(world.sf)
      if (world.sf.num.before == world.sf.num.after){
        cat(paste0("kriging: Cannot find the country: ", cty[q],"\n"))
      }
    }
    #add back subresgions to mask
    addback.map.mask <- c()
    for (k in 1:length(cty)){
      if (grepl("japan", cty[k], ignore.case = T)){
        if (length(exclude.region) != 0){
          exclude.region <- unlist(strsplit(exclude.region, ','))
          addback.map.mask <- st_as_sf(maps::map(cty[k], exclude.region, plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) %>% st_geometry() %>% st_cast('POLYGON') %>% st_union()
        }
      }
    }
    
    colnames(subset.out.table.kriged)[3] <- "Longitude"
    colnames(subset.out.table.kriged)[4] <- "Latitude"
    if (grepl("light", krig.map.theme)){
      #light theme
      out <- gsub("txt$", paste0(target, "_vs_",pops[n],".kriging.light_theme.tiff"), info)
      #out <- gsub("txt$", paste0(trait.name, ".kriging.light_theme.tiff"), file)
      #out <- gsub("\\.+",".", out)
      theme.country.line.color <- "black"
      theme.country.fill <- "white"
      theme.background <- "lightgray"
      theme.color.low <- "gold2"
      theme.color.high <- "blue"
    } else if (grepl("dark", krig.map.theme)){
      #dark theme
      out <- gsub("txt$", paste0(target, "_vs_",pops[n],".kriging.dark_theme.tiff"), info)
      #out <- gsub("txt$", paste0(trait.name, ".kriging.dark_theme.tiff"), file)
      #out <- gsub("\\.+",".", out)
      theme.country.line.color <- "grey50"
      theme.country.fill <- "black"
      theme.background <- "grey85"
      theme.color.low <- "gold2"
      theme.color.high <- "purple1"
    } else {
      #rainbow theme
      out <- gsub("txt$", paste0(target, "_vs_",pops[n],".kriging.rainbow_theme.tiff"), info)
      #out <- gsub("txt$", paste0(trait.name, ".kriging.rainbow_theme.tiff"), file)
      #out <- gsub("\\.+",".", out)
      theme.country.line.color <- "grey50"
      theme.country.fill <- "grey30"
      theme.background <- "grey95"
    }
    out.plot <-
      ggplot() +
      labs(title = paste0(target, " vs. ",pops[n])) +
      geom_point(data = as.data.frame(subset.out.table.kriged), aes(Longitude, Latitude, color=var1.pred), size = kg/length(grd), shape = 15) +
      geom_sf(data = region.map.mask, 
              fill = theme.background, #background color
              color = NA) +
      geom_sf(data = world.sf,
              fill = theme.country.fill, #country color
              color = NA)
    if (length(addback.map.mask) != 0){
      out.plot <- out.plot +
        geom_sf(data = addback.map.mask,
                fill = theme.country.fill, #country color
                color = NA)        
    }
    out.plot <- out.plot +
      theme(text = element_text(size = text.size),
            panel.background = element_rect(fill = theme.background, color = theme.background), #background color
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            #panel.border = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
            axis.text = element_text(color = "black", size = text.size-2),
            axis.title = element_text(color = "black", size = text.size),
            axis.ticks = element_line(colour = "black", size = 0.2),
            #legend.title = element_text(size = text.size-2),
            #legend.title = element_blank(),
            legend.text = element_text(size = text.size),
            legend.key = element_rect(fill = NA, color = NA),
            legend.key.height= unit(fig_height*0.1, 'cm'),
            legend.key.width= unit(fig_width*0.05, 'cm'),
            plot.title = element_text(hjust = 0.5))
    if (grepl("light|dark", krig.map.theme)){
      world.sf.line <- st_as_sf(maps::map("world2", plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) #load world map
      out.plot <- out.plot +
         #geom_sf(data = world.sf.line,
          #       fill = NA,
          #       color = theme.country.line.color,
          #       size = 0.1) +
        scale_color_gradient(low=theme.color.low,
                             high=theme.color.high,
                             name = "Genetic\ndistance")
    } else {
      if(!require("viridisLite")) install.packages("viridisLite")
      library(viridisLite)
      world.sf.line <- st_as_sf(maps::map("world2", plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) #load world map
      out.plot <- out.plot +
        #geom_sf(data = world.sf.line,
        #        fill = NA,
        #        color = theme.country.line.color,
        #        size = 0.2) +
        scale_color_viridis_c(name = "Genetic\ndistance",
                              option = "H",
                              direction = -1)
    }
    out.plot <- out.plot +
      coord_sf(xlim = c(left,right), ylim = c(bottom,top))
    if (length(leg.pos) > 0){
      if (leg.pos == "bottomright"){
        out.plot <- out.plot +
          theme(
            legend.justification = c(1, 0), 
            legend.position = c(1, 0)
          )
      }
      if (leg.pos == "topright"){
        out.plot <- out.plot +
          theme(
            legend.justification = c(1, 1), 
            legend.position = c(1, 1)
          )
      }
      if (leg.pos == "bottomleft"){
        out.plot <- out.plot +
          theme(
            legend.justification = c(0, 0), 
            legend.position = c(0, 0)
          )
      }
      if (leg.pos == "topleft"){
        out.plot <- out.plot +
          theme(
            legend.justification = c(0, 1), 
            legend.position = c(0, 1)
          )
      }
    }
    #pdf(out, width=fig_width, height=fig_height)
    tiff(out, units="cm", width = fig_width, height = fig_height, res = 600)
    print(out.plot)
    dev.off()
  }
}
cat("Done.\n")



