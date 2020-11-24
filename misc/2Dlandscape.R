#!/usr/bin/Rscript --vanilla

# list of required libraries
list.of.packages <- c("ggplot2","akima","optparse")

# check if required libraries are installed and quit otherwise
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {
  print("The following R packages are missing to properly execute this script:")
  for (p in new.packages) {
    print(p)
  }
  print("Please install the above packages and re-run the script")
  quit("no", status = 1)
}

# load require libraries
for (package in list.of.packages) {
  library(package, character.only = T)
}

# command line argument specifications
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="RNA2Dfold input file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="2Dplot.pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p","--partfunc"), type="logical", action="store_true", default=FALSE,
              help="Plot Ensemble Energy instead of MFE, i.e. use 7th column instead of 6th [default= %default]"),
  make_option(c("--maxX"), type="integer", default=-1,
              help="maximum X value"),
  make_option(c("--maxY"), type="integer", default=-1,
              help="maximum Y value"),
  make_option(c("--minZ"), type="integer", default=10000,
              help="minimum Z value, i.e. minimum energy"),
  make_option(c("--maxZ"), type="integer", default=0,
              help="maximum Z value, i.e. maximum energy"),
  make_option(c("--binX"), type="integer", default=100,
              help="number of bins in X dimension [default= %default]"),
  make_option(c("--binY"), type="integer", default=100,
              help="number of bins in Y dimension [default= %default]"),
  make_option(c("--labX"), type="character", default="Distance to 1st Reference",
              help="Label for X-axis, i.e. 1st reference structure", metavar="character"),
  make_option(c("--labY"), type="character", default="Distance to 2nd Reference",
              help="Label for Y-axis, i.e. 2nd reference structure", metavar="character"),
  make_option(c("--Xres"), type="integer", default=10,
              help="Resolution for X values"),
  make_option(c("--Yres"), type="integer", default=10,
              help="Resolution for Y values"),
  make_option(c("--Zres"), type="integer", default=4,
              help="Resolution for Z value iso-lines, i.e. energies")
);

opt_parser = OptionParser(usage           = "usage: %prog [options]",
                          description     = "An R script to create colored 2D landscape heatmaps from RNA2Dfold output.",
                          add_help_option = TRUE,
                          option_list     = option_list);

opt = parse_args(opt_parser);

# quit if input file is missing
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Specify at lease the input file", call.=FALSE)
}

# by default, we skip the first 6 lines of RNA2Dfold output
lines_to_skip <- 6

# read first line of input file to see whether we need to skip
# the FASTA header as well
f = file(opt$file,open="r")
line = readLines(f, n = 1)
if (substr(line, 1, 1) == '>') {
  lines_to_skip <- 7
}
close(f)

# read input data
dat <- read.table(opt$file, head=F,skip=lines_to_skip)

colnames(dat)[colnames(dat)=="V1"] <- "x"
colnames(dat)[colnames(dat)=="V2"] <- "y"
if (opt$partfunc) {
  colnames(dat)[colnames(dat)=="V7"] <- "z"
} else {
  colnames(dat)[colnames(dat)=="V6"] <- "z"
}

# remove Inf values
dat <- dat[!(dat$z == "Inf"),]

# compute missing values if not specified through command line options

if (opt$maxX == -1) {
  opt$maxX <- opt$Xres * ceiling(max(dat$x) / opt$Xres)
}

if (opt$maxY == -1) {
  opt$maxY <- opt$Yres * ceiling(max(dat$y) / opt$Yres)
}

if (opt$minZ == 10000) {
  opt$minZ <- opt$Zres * floor(min(dat$z) / opt$Zres)
}

# start preparing data
if (max(dat$z) > opt$maxZ) {
  dat[(dat$z > opt$maxZ),]$z = opt$maxZ
}

if (min(dat$z) < opt$minZ) {
  dat[(dat$z < opt$minZ),]$z = opt$minZ
}

di <- interp(dat$x, dat$y, dat$z,
             xo=seq(0, opt$maxX, length=opt$binX),
             yo=seq(0, opt$maxY, length=opt$binY))

dat_interp <- data.frame(expand.grid(x=di$x, y=di$y), z=c(di$z))


contour_breaks_strong = seq(opt$minZ, opt$maxZ, opt$Zres)
contour_breaks_weak   = seq(opt$minZ - ( opt$Zres / 2.), opt$maxZ - ( opt$Zres / 2.), opt$Zres)

p <- ggplot(dat_interp, aes(x=x, y=y, z=z, fill=z))
p <- p +  geom_raster(na.rm=T, alpha=1, interpolate=F)
p <- p +  stat_contour(aes(z=z), colour="black", size=0.2, breaks=contour_breaks_strong)
p <- p +  stat_contour(aes(z=z), colour="black", size=0.05, linetype = "dotted", breaks=contour_breaks_weak)
p <- p +  scale_fill_distiller( type="div",
                        palette = 7,
                        limits  = c(opt$minZ, opt$maxZ),
                        breaks = seq(opt$minZ, opt$maxZ, 10),
                        na.value="transparent")
p <- p +  scale_x_continuous( limits=c(0, opt$maxX),
                              expand = c(0, 0),
                              breaks=seq(opt$Xres, opt$maxX - opt$Xres, opt$Xres))
p <- p +  scale_y_continuous( limits=c(0, opt$maxY),
                              expand = c(0, 0),
                              breaks=seq(opt$Yres, opt$maxY - opt$Yres, opt$Yres))
p <- p +  xlab(opt$labX)
p <- p +  ylab(opt$labY)
p <- p +  theme_bw()
p <- p +  theme(plot.background   = element_blank(),
                panel.border      = element_blank(),
                panel.grid.major  = element_blank(),
                panel.grid.minor  = element_blank(),
                axis.line = element_line(colour = "black"))
p <- p +  theme(axis.title.x  = element_text(family="Helvetica", colour="#000000"),
                axis.text.x   = element_text(family="Helvetica", colour="#111111"),
                axis.title.y  = element_text(family="Helvetica", colour="#000000"),
                axis.text.y   = element_text(family="Helvetica", colour="#111111"),
                legend.position = "top"
               )
p <- p +  guides(fill=guide_colorbar(title="Energy (kcal/mol)", title.position="top", title.vjust=0.8, barwidth = 20, barheight = 1))


ggsave(filename=opt$out, plot=p, width=5, height=6)
