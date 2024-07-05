#install.packages("vctrs", repos = "https://packagemanager.rstudio.com/cran/latest")
#install.packages("scatterpie", repos = "https://packagemanager.rstudio.com/cran/latest")
#install.packages("extrafont", repos = "https://packagemanager.rstudio.com/cran/latest")
#usage: Rscript Global_map_ggplot.R data_file region_border figure_size dot_size
#region_border, figure_size and dot_size are optional
#region border format: left,right,bottom,top (eg. 70,150,15,57), if do kriging, you must contain entire country in the map
#google format will be transformed to: (eg. 70,15,150,57)
#Stadia API_key: XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX
#figure size format: width,height (eg. 10,7)
options(warn=-1)
s.help <- function(){
  cat("\nThis scirpt is writtern by Ben Chien. Jul. 2022
Usage: Rscript global_map_ggplot.R -f TABLE_FILE [-d] [-p] [-vc INT] [-cc INT] [-r MAP_BORDER] [-t API_KEY] [-fs WIDTH_HEIGHT] [-tm THEME] [-lp POSITION] [-k NUMBER] [-ktm THEME] [-cty NAME] [--kmodel MODEL] [--kpsill NUMBER] [--knugget NUMBER] [--krange NUMBER] [--kregion MAP_BORDER] [-h]\n\n
-t/--terrain: using terrain map.
-p/--pie: pie chart.
-f/--file: table file. Table seperates by tab with header line.
-d/--discrete: implying discrete type of dataset. default: false
-cc/--coord_column: implying column index that contains coordinates (longitude). Latitude should be next to longitude column. default: 2
-vc/--value_column: implying column index that contains values. If defined by 2 columns, the second one will present as shapes. eg. 4,6. default: 4
-r/--region: map coordinate border (left,right,bottom,top). default: -180,180,-90,90
-fs/--figure_size: figure size (width,height). default: 10,6
-tm/--theme: map theme. Possible values: light, dark. default: light
-lp/--legend_position: legend position. Possible values: topright, bottomright, topleft, bottomleft, none. default: system default
-k/--kriging_res: Kriging resolution in kilometer. default: 10. Be careful: smaller number may take more time.
  **You can increase the resolution at the first time to have a quick look of the fitting model. 
-cty/--country: country name that you are looking at. Multiple regions can be seperated by comma.
-ktm/--ktheme: map theme. Possible values: rainbow, light, dark. default: rainbow
--kmodel: model used for Krigging. default: Sph
--kpsill: psill value used for Krigging. default: interval formula
--knugget: nugget value used for Krigging. default: interval formula
--krange: range value used for Krigging. default: interval formula
--kregion: map coordinate border (left,right,bottom,top) for kriging plot. default: -180,180,-90,90
  **If you got incorrect output, increasing border range of --kregion might solve the problem.
-h/--help: help.\n\n")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  s.help()
  quit()
}
terrain <- c()
pie <- 0 #pie chart
file <- c()
region <- c()
fig_size <- c()
leg.pos <- c()
leg.pos.values <- c("bottomright", "bottomleft", "topright", "topleft")
dis <- 0 #discrete type of data, true or false
col <- 4 #column index that contains values for plotting
col.2 <- c()
coord.col <- 2 #column index that contains coordinates (long) for plotting
map.theme <- "light"
kg <- c()
krig.map.theme <- "rainbow"
cty <- c()
krig.range.in <- c()
krig.nugget.in <- c()
krig.psill.in <- c()
krig.model <- c()
krig.region <- c()
for (i in 1:length(args)){
  if (args[i] == '-f' || args[i] == '--file'){
    file <- args[i+1]
    if (!file.exists(file)){
      s.help()
      cat("-f: file does not exist.\n")
      quit()
    }
  }
  if (args[i] == '-h' || args[i] == '--help'){
    s.help()
    quit()
  }
  if (args[i] == '-t' || args[i] == '--terrain'){
    terrain <- args[i+1]
  }
  if (args[i] == '-p' || args[i] == '--pie'){
    pie <- 1
  }
  if (args[i] == '-r' || args[i] == '--region'){
    region <- args[i+1]
    #128,150,30,46 for japan [22*16=352] zoom = 7, 22+16=38 --> 20*20
    #70,150,15,57 for V. angularis [80*42=3360] zoom = 5, 80+42=122 --> 60*60
  }
  if (args[i] == '--kregion'){
    krig.region <- args[i+1]
    #128,150,30,46 for japan
    #70,150,15,57 for V. angularis
  }
  if (args[i] == '-fs' || args[i] == '--figure_size'){
    fig_size <- args[i+1]
    #4.2,3 for japan
    #10,7 for V. angularis
  }
  if (args[i] == '-d' || args[i] == '--discrete'){
    dis <- 1
  }
  if (args[i] == '-vc' || args[i] == '--value_column'){
    col <- args[i+1]
    if (grepl(',', col)){
      col <- unlist(strsplit(col, ','))
      if (grepl("[^0-9]", col[2])){
        s.help()
        cat("-c: value is wrong.\n")
        quit()
      }
      col.2 <- as.integer(col[2])
    }
    if (grepl("[^0-9]", col[1])){
      s.help()
      cat("-c: value is wrong.\n")
      quit()
    }
    col <- as.integer(col[1])
  }
  if (args[i] == '-cc' || args[i] == '--coord_column'){
    coord.col <- args[i+1]
    if (grepl("[^0-9]", coord.col)){
      s.help()
      cat("-cc: value is wrong.\n")
      quit()
    }
    coord.col <- as.integer(coord.col)
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
  if (args[i] == '-k' || args[i] == '--kriging_res'){
    if (!grepl("^-", args[i+1]) && i < length(args)){
      kg <- as.numeric(args[i+1])
      #you can try 10 first
      kg <- kg * 1000
    } else {
      kg <- 10000
    }
  }
  if (args[i] == '-ktm' || args[i] == '--ktheme'){
    if (grepl("light|dark|rainbow", args[i+1])){
      krig.map.theme <- args[i+1]
    } else {
      cat("-ktm: cannot recognize the theme used. Default setting used.\n")
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
}

if (length(file) < 1){
  s.help()
  cat("-f: data file is not provided, only draw a blank map.\n")
  length(kg) = 0
  pie = 0
}

if (length(terrain) > 0){
  library(ggmap)
  register_stadiamaps(terrain)
  map.theme <- "light"
}

if (pie == 1){
  library(scatterpie)
  dis <- 1
}

if (length(kg) > 0){
  if(!require("gstat")) install.packages("gstat")
  if(!require("geosphere")) install.packages("geosphere")
  if(!require("sp")) install.packages("sp")
  library(gstat)
  library(geosphere)
  library(sp)
  if (length(krig.model) == 0){
    krig.model <- "Sph"
  }
  if (length(krig.range.in) != 0 || length(krig.nugget.in) != 0 || length(krig.psill.in) != 0){
    if (length(kg) == 0){
      cat("-k: Kriging option is not set. --krange, --kpsill or --knugget will not be used.\n")
    }
  }
}

if(!require("ggplot2")) install.packages("ggplot2")
if(!require("sf")) install.packages("sf")
if(!require("maps")) install.packages("maps")
library(ggplot2)
library(sf)
library(maps)
library(extrafont)
loadfonts(device = "win", quiet = TRUE)

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
  if (length(terrain) > 0){
    t_region <- c(left, bottom, right, top)
    #128,150,30,46 for japan [22*16=352] zoom = 7, 22+16=38 --> 20*20
    #70,150,15,57 for V. angularis [80*42=3360] zoom = 5, 80+42=122 --> 60*60
    t_area <- sqrt((right-left)*(top-bottom))
    zoom <- 5
    if (t_area <= 2){
      zoom <- 10
    } else if (t_area <= 5 && t_area > 2){
      zoom <- 9
    } else if (t_area <= 10 && t_area > 5){
      zoom <- 8
    } else if (t_area <= 20 && t_area > 10){
      zoom <- 7
    } else if (t_area <= 40 && t_area > 20){
      zoom <- 6
    } else if (t_area <= 70 && t_area > 40){
      zoom <- 5
    } else if (t_area <= 110 && t_area > 70){
      zoom <- 4
    } else if (t_area <= 160 && t_area > 110){
      zoom <- 3
    } else if (t_area <= 220 && t_area > 160){
      zoom <- 2
    } else if (t_area > 220){
      zoom <- 1
    }
    t_map <- get_stadiamap(t_region, zoom = zoom, maptype = "stamen_terrain_background")
  } 
}
if (bottom == -90){
  bottom = -89.9999
}
if (top == 90){
  top = -89.9999
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

#load world map
world.sf <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>% st_transform(crs = 4326) #load world map
#read data
my.data <- c()
if (length(file) != 0){
  my.data <- read.table(file, header = T, as.is = T, sep = "\t")
  my.data <- my.data[complete.cases(my.data[,col]),] #remove NA data
}
if (dis == 1 && length(file) != 0){
  #colors <- c("royalblue3", "tan3", "gold2", "purple4", "navy", "red4", "dodgerblue4", "tan1", "orange3", "deeppink2", "purple1")
  colors <- c("gold","#ed3325","#4b8bcb","#c6b1d4","#255271", "#f7931e", "#ed7f6d", "#7758a5", "#f9bbb9", "#c6b1d4", "#921a1d", "#90a7b7")
  my.data[,col] <- as.character(my.data[,col])
  if (length(col.2) != 0){
    my.data <- my.data[complete.cases(my.data[,col.2]),] #remove NA data
    my.data[,col.2] <- as.character(my.data[,col.2])
  }
  dis.cats <- length(unique(my.data[,col]))
  if (pie == 1){
    dis.cats <- length(unique(my.data[1,])) - 4 #do not count: group(ID), long, lat, total_size
  }
  if (dis.cats > 11){
    if(!require("RColorBrewer")) install.packages("RColorBrewer")
    library(RColorBrewer)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    colors <- getPalette(dis.cats)
  }
}
#set label markers
#breaks <- c(min(my.data[,4]),(max(my.data[,4])-min(my.data[,4]))/2 , max(my.data[,4]))
#plot
legend.name <- c()
trait.name <- "blank"
output_name <- paste0(getwd(), "/map.txt")
if (length(file) != 0){
  legend.name <- colnames(my.data)[col]
  legend.name <- gsub("\\.+"," ", legend.name)
  legend.name <- gsub(" $","", legend.name)
  trait.name <- colnames(my.data)[col]
  output_name <- file
}
if (grepl("light", map.theme)){
  #light theme
  out <- gsub("txt$", paste0(trait.name, ".light_theme.pdf"), output_name)
  out <- gsub("\\.+",".", out)
  theme.country.line.color <- "black"
  theme.country.fill <- "white"
  theme.background <- "lightgray"
  theme.color.low <- "gold2"
  theme.color.high <- "blue"
} else {
  #dark theme
  out <- gsub("txt$", paste0(trait.name, ".dark_theme.pdf"), output_name)
  out <- gsub("\\.+",".", out)
  theme.country.line.color <- "grey50"
  theme.country.fill <- "black"
  theme.background <- "grey85"
  theme.color.low <- "gold2"
  theme.color.high <- "purple1"
}
if (length(terrain) > 0){
  out <- gsub("txt$", paste0(trait.name, ".terrain_theme.pdf"), output_name)
  out <- gsub("\\.+",".", out)
  theme.background <- "white"
  theme.country.fill <- NA
  gg.out <- ggmap(t_map, darken = c(0.5, "white"), extent = "device")# + coord_fixed(ratio = 1)
} else {
  gg.out <- ggplot() + 
    geom_sf(data = world.sf,
            color = theme.country.line.color, 
            fill = theme.country.fill, 
            size = 0.1) +
    xlim(left,right) + 
    ylim(bottom,top)
}
if (pie == 1){
  data_start <- col
  data_end <- ncol(my.data) - 1 #don't use total_size
  if (coord.col > col){
    data_end <- coord.col - 2
  }
  my.data[,data_start:(data_end+1)] <-sapply(my.data[,data_start:(data_end+1)],function(x) as.numeric(as.character(x)))
  colnames(my.data)[1] <- "group"
  colnames(my.data)[coord.col] <- "long"
  colnames(my.data)[(coord.col+1)] <- "lat"
  colnames(my.data)[(data_end+1)] <- "total"
  gg.out <- gg.out + 
    geom_scatterpie(data = my.data,
                    aes(x = long, 
                        y = lat,
                        group = group,
                        r = log10(total)*2,
                    ),
                    cols = colnames(my.data)[data_start:data_end],
                    alpha = 0.9)
  #breaks <- c(log10(10)*2,log10(100)*2,log10(500)*2)
} else {
  if (length(file) != 0 && length(col.2) == 0){
    gg.out <- gg.out + geom_point(data = my.data, 
                                  aes(my.data[,coord.col], 
                                      my.data[,coord.col+1], 
                                      color = my.data[,col]), 
                                  alpha =0.85, 
                                  shape = 16, 
                                  size = round(0.2*min(c(fig_width, fig_height), na.rm = T), digits = 1))
  } else if (length(file) != 0) {
    gg.out <- gg.out + geom_point(data = my.data, 
                                  aes(my.data[,coord.col], 
                                      my.data[,coord.col+1], 
                                      color = my.data[,col],
                                      shape = my.data[,col.2]),
                                  alpha =0.85,
                                  size = round(0.4*min(c(fig_width, fig_height), na.rm = T), digits = 1))  
  }
}
gg.out <- gg.out + theme(text = element_text(family="sans", 
                                             size = round(3.5*min(c(fig_width, fig_height), na.rm = T),digits = 0)
                                             #size = 9/.pt
                                             ),
        panel.background = element_rect(fill = theme.background, color = theme.background),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        legend.text = element_text(size = round(3.5*min(c(fig_width, fig_height)-3, na.rm = T), digits = 0)),
        legend.title = element_text(face = "bold", size = round(3.5*min(c(fig_width, fig_height)-2, na.rm = T), digits = 0)),
        legend.key = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = alpha(theme.background, 0.4), color = alpha(theme.background, 0.4))
  ) + 
  xlab("Longitude") + 
  ylab("Latitude")
if (length(file) != 0 && dis == 0){
  gg.out <- gg.out +
    scale_color_gradient(low=theme.color.low, 
                         high=theme.color.high, 
                         name = legend.name, 
                         #breaks = breaks, 
                         #labels = format(breaks)
    )
} else {
  if (pie == 1){
    gg.out <- gg.out +
      scale_fill_manual(values = colors) +
      theme(panel.border = element_rect(fill = NA, color = 'black', size = 1),
            legend.background = element_rect(fill = alpha(theme.background, 0.4), color = NA)
            )
  } else if (length(file) != 0) {
    gg.out <- gg.out +
      scale_color_manual(values = colors,
                         name = legend.name, 
                         #breaks = breaks, 
                         #labels = format(breaks)
      )
  }
  if (length(file) != 0 && length(col.2) != 0){
    trait.name2 <- gsub("\\.+"," ", colnames(my.data)[col.2])
    trait.name2 <- gsub(" $","", trait.name2)
    gg.out <- gg.out +
      guides(shape = guide_legend(title = trait.name2))
  }
}
if (length(leg.pos) > 0){
  if (leg.pos == "bottomright"){
    gg.out <- gg.out +
      theme(
        legend.justification = c(1, 0), 
        legend.position = c(1, 0)
      )
  }
  if (leg.pos == "none"){
    gg.out <- gg.out +
      theme(
        legend.position = 'none'
      )
  }
  if (leg.pos == "topright"){
    gg.out <- gg.out +
      theme(
        legend.justification = c(1, 1), 
        legend.position = c(1, 1)
      )
  }
  if (leg.pos == "bottomleft"){
    gg.out <- gg.out +
      theme(
        legend.justification = c(0, 0), 
        legend.position = c(0, 0)
      )
  }
  if (leg.pos == "topleft"){
    gg.out <- gg.out +
      theme(
        legend.justification = c(0, 1), 
        legend.position = c(0, 1)
      )
  }
}
pdf(out, width=fig_width, height=fig_height)
print(gg.out)
dev.off()

#kriging plot
if (length(kg) > 0 && dis == 1){
  cat("-kg: Discrete type of dataset cannot do kriging.\n")
}
if (length(kg) > 0 && (length(terrain) > 0 || pie == 1)){
  cat("-kg: It is incompatible to use -kg and -t/-p at the same time.\n")
}
if (length(kg) > 0 && dis == 0){
  if (length(krig.left) == 0){
    krig.left <- left
    krig.right <- right
    krig.bottom <- bottom
    krig.top <- top
  }
  if (col > coord.col){
    data.idx <- col - 2
  } else {
    data.idx <- col
  }
  if (length(cty) != 0){
    cty <- unlist(strsplit(cty, ','))
  } else {
    cty <- '.'
  }
  #my.data.ori <- my.data
  colnames(my.data)[coord.col] <- "long"
  colnames(my.data)[coord.col+1] <- "lat"
  coordinates(my.data) <- ~long+lat
  #remove duplicate coord first
  my.data <- my.data[which(!duplicated(my.data@coords)),]
  #transform coordinate system to meter unit
  proj4string(my.data) <- CRS("+init=epsg:4326")
  CRS.new <- CRS("+init=epsg:3857") #meter coord system
  my.data.meter <- spTransform(my.data, CRS.new)
  #transform data to vgm variables
  my.data.vgm <- variogram(my.data.meter@data[[data.idx]]~1, my.data.meter)
  #my.data.vgm.test <- plot(variogram(my.data.meter@data[[data.idx]]~1, my.data.meter, alpha = c(0, 45, 90, 135)))
  #out <- gsub("txt$", paste0(trait.name, ".k_model.test.pdf"), file)
  #out <- gsub("\\.+",".", out)
  #pdf(out, width=4, height=3)
  #print(my.data.vgm.test)
  #dev.off()
  #remove any duplicate coordinate
  #my.data.meter <- my.data.meter[which(!duplicated(my.data.meter@coords)),]
  #set up Kriging parameters
  if (length(krig.range.in) == 0){
    krig.range <- round(max(my.data.vgm$dist)-max(my.data.vgm$dist)/100)
  } else {
    krig.range <- krig.range.in
  }
  krig.sill.digit <- ceiling(-log10(max(my.data.vgm$gamma)))
  krig.sill <- quantile(my.data.vgm$gamma, 0.9)
  if (length(krig.nugget.in) == 0){
    krig.nugget <- as.numeric(quantile(my.data.vgm$gamma[1:5], 0.1))-((krig.sill-as.numeric(quantile(my.data.vgm$gamma[1:5], 0.1)))/length(my.data.vgm$gamma))
  } else {
    krig.nugget <- krig.nugget.in
  }
  if (krig.nugget < 0){
    krig.nugget = 0
  }
  if (length(krig.psill.in) == 0){
    if (krig.sill - krig.nugget > 0){
      krig.psill <- krig.sill - krig.nugget
    } else {
      cat("--knugget value is too big. Reseted by the interval formula.\n")
      krig.nugget <- as.numeric(quantile(my.data.vgm$gamma[1:5], 0.1))-((krig.sill-as.numeric(quantile(my.data.vgm$gamma[1:5], 0.1)))/length(my.data.vgm$gamma))
      krig.psill <- krig.sill - krig.nugget
    }
  } else {
    krig.psill <- krig.psill.in
  }
  #apply model (Exp, Sph, Lin or Gau)
  model.krig <- vgm(krig.psill, model=krig.model, nugget=krig.nugget, range=krig.range)
  #fit model
  fit.variog <- fit.variogram(my.data.vgm, model.krig, fit.kappa = T, warn.if.neg = F)
  #plot fitness of the model
  fitness <- plot(my.data.vgm, fit.variog, main = legend.name)
  out <- gsub("txt$", paste0(trait.name, ".k_model.fitness.pdf"), file)
  out <- gsub("\\.+",".", out)
  pdf(out, width=4, height=3)
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
  x.range <- as.integer(c(x.range.min,x.range.max))
  y.range <- as.integer(c(y.range.min,y.range.max))
  grd <- expand.grid(long=seq(from=x.range[1], to=x.range[2], by=kg), lat=seq(from=y.range[1], to=y.range[2], by=kg))
  coordinates(grd) <- ~long+lat
  proj4string(grd) <- CRS("+init=epsg:3857")
  gridded(grd) <- TRUE
  #do Kriging modelling
  cat("Applying Kriging model to the map...\nIt may take some time...\n")
  my.data.kriged <- krige(my.data.meter@data[[data.idx]]~1, my.data.meter, grd, model = model.krig)
  proj4string(my.data.kriged) <- CRS("+init=epsg:3857")
  my.data.kriged <- spTransform(my.data.kriged, CRS("+init=epsg:4326"))
  my.data.kriged@data[["var1.pred"]] <- my.data.kriged@data[["var1.pred"]]
  my.data.kriged@data[["var1.var"]] <- my.data.kriged@data[["var1.var"]]
  my.data.kriged <- as.data.frame(my.data.kriged)
  #do masking of the rest region except target region
  region.map.mask <- st_as_sf(maps::map("world", cty, plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) %>% st_geometry() %>% st_cast('POLYGON') %>% st_union()
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
  colnames(my.data.kriged)[3] <- "Longitude"
  colnames(my.data.kriged)[4] <- "Latitude"
  if (grepl("light", krig.map.theme)){
    #light theme
    out <- gsub("txt$", paste0(trait.name, ".kriging.light_theme.pdf"), file)
    out <- gsub("\\.+",".", out)
    theme.country.line.color <- "black"
    theme.country.fill <- "white"
    theme.background <- "lightgray"
    theme.color.low <- "gold2"
    theme.color.high <- "blue"
  } else if (grepl("dark", krig.map.theme)){
    #dark theme
    out <- gsub("txt$", paste0(trait.name, ".kriging.dark_theme.pdf"), file)
    out <- gsub("\\.+",".", out)
    theme.country.line.color <- "grey50"
    theme.country.fill <- "black"
    theme.background <- "grey85"
    theme.color.low <- "gold2"
    theme.color.high <- "purple1"
  } else {
    #rainbow theme
    out <- gsub("txt$", paste0(trait.name, ".kriging.rainbow_theme.pdf"), file)
    out <- gsub("\\.+",".", out)
    theme.country.line.color <- NA
    theme.country.fill <- "#30123BFF"
    theme.background <- "grey95"
  }
  out.plot <-
    ggplot() +
    geom_point(data = as.data.frame(my.data.kriged),aes(Longitude, Latitude, color=var1.pred), size = kg/length(grd), shape = 15) +
    geom_sf(data = region.map.mask, 
            fill = theme.background, #background color
            color = NA) +
    geom_sf(data = world.sf,
            fill = theme.country.fill, #country color
            color = NA) +
    theme(text = element_text(size = round(3.5*min(c(fig_width, fig_height), na.rm = T), digits = 0)),
          panel.background = element_rect(fill = theme.background, color = theme.background), #background color
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(face = "bold", size = round(3.5*min(c(fig_width, fig_height)-2, na.rm = T), digits = 0)),
          legend.text = element_text(size = round(3.5*min(c(fig_width, fig_height)-3, na.rm = T), digits = 0)),
          legend.key = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  if (grepl("light|dark", krig.map.theme)){
    world.sf.line <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>% st_set_crs(4326) #load world map
    out.plot <- out.plot +
      geom_sf(data = world.sf.line,
              fill = NA,
              color = theme.country.line.color,
              size = 0.1) +
      scale_color_gradient(low=theme.color.low,
                          high=theme.color.high,
                          name = legend.name)
  } else {
    if(!require("viridisLite")) install.packages("viridisLite")
    library(viridisLite)
    out.plot <- out.plot +
      scale_color_viridis_c(name = legend.name,
                            option = "H",
                            direction = 1)
  }
  out.plot <- out.plot +
    coord_sf(xlim = c(left,right), ylim = c(bottom,top))
  if (length(leg.pos) > 0){
    if (leg.pos == "none"){
      out.plot <- out.plot +
        theme(
          legend.position = 'none'
        )
    }
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
  pdf(out, width=fig_width, height=fig_height)
  print(out.plot)
  dev.off()
}

