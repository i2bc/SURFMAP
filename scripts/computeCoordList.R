rm(list=ls())


###################  LIBRARIES   ##########################

library(optparse)
library(data.table)

###########################################################



################## FUNCTIONS ##########################

# This function returns a data frame representing the subdivision of the energy 
# map into squares of equal area. The size of the squares is determined by the
# integers quadrillage1 and quadrillage 2. For each square several values are
# calculated, such as : 
# - mean Energy of the points within the square
# - number of points within the square
# - theta/phi coordinates of the square
# - normalized (compared to lowest energy in the data frame) mean energy

energy_matrix <- function(data, nb_file) {
  # Determine number of squares ie number of element of the vector which will 
  # serve for the energy matrix. 
  # It is necessary that : 360 % quadrillage1 = 0 & 180 % quadrillage2 = 0
  data[,6] <- seq(1, nrow(data), by = 1)
  corres = data[,c(4,5,6)]
  data <- data[,-c(4,5)]
  seqcol = seq(-180,180,by=width)
  seqrow = seq(-90,90,by=height)
  
  # First using sapply to gather all coordinate belonging to a same pixel. Result is a list of list of dataframe/
  # First level list contain list of dataframe of elem belonging to a same column (for ex. all pixel of abscissa coords belonging to (80,90)). These element are grouped into dataframes.
  datasp1 = split(data,findInterval(data[,1], seqcol))
  nabs <- names(datasp1)
  datasp1 = lapply(nabs, function(x) {datasp1[[x]][,1] = as.integer(x)*width; return(datasp1[[x]])})
  datasp2 = lapply(datasp1, function(x) split(x, findInterval(x[,2], seqrow)))
  datasp3 = sapply(datasp2, function(x) {nord = names(x); 
                   lapply(nord, function(y) {x[[y]][,2] = as.integer(y)*height; return(x[[y]])})})
  datasp4 = unlist(datasp3, recursive = FALSE)
  dataf = do.call(rbind,datasp4)

  colnames(dataf)[1] = "abscissa"
  colnames(dataf)[2] = "ordinate"
  #colnames(dataf)[1] = "phi"
  #colnames(dataf)[2] = "theta"
  colnames(dataf)[3] = "score"
  colnames(dataf)[4] = "index_sol"

  colnames(corres)[1] = "resnb"
  colnames(corres)[2] = "restype"
  colnames(corres)[3] = "index_sol"

  dataf = merge(dataf, corres, "index_sol")
  datafinal = dataf[,c(2,3,4,5,6,1)]
  datafinal = datafinal[order(datafinal[,1], datafinal[,2]),]

  return(datafinal)
}

#####################################################



#####################################################
#                      MAIN
#####################################################

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NA, 
              help="path to files", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NA, 
              help="(phi, theta, val) file", metavar="character"),
  make_option(c("-s", "--cellsize"), type="character", default=5, 
              help="dimension of a grid cell. max is 180. please chose a multiple of 180. (default = 5)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (!is.na(opt$file)) {
  path = dirname(opt$file)
  files = c(basename(opt$file))
} else if (!is.na(opt$path)) {
  path = opt$path
  files = list.files(path, pattern = "\\partlist.out$")
} else {
  cat("\nError: no file nor directory in input. Exiting now.\n")
  stop()
}

# Size of the grid.
if (!is.na(opt$cellsize)) {
  width <<- as.double(opt$cellsize)
} else {
  width <<- 5
}
height <<- width

if (180%%width != 0) {
  cat("\nError: Grid cell is not a multiple of 180. Exiting now.\n")
  stop()
}

setwd(path)
dir.create("coord_lists", showWarnings = FALSE)

for (file in (1:length(files))) {
  name_prefix = gsub("_partlist.out", "", files[file]) 
  #cat(name_prefix, "\n")
  data = read.table(files[file], fill = TRUE, header = FALSE)
  
  data[,1] = (data[,1])%%(2*pi) # Now phi values are included in [0, 2*pi]
  data[,1] = data[,1]-pi # Now phi values are included in [-pi, pi]

  # sinusoidal projection
  # If theta belongs to [0, pi], phiproj=phi*sin(theta),
  # If theta belongs to [-pi/2, pi/2], phiproj=phi*cos(theta)
  data[,1] = data[,1]*sin(data[,2]) 
  data[,1] = data[,1]*180/pi
  data[,2] = data[,2]*180/pi
  data[,2] = 90 - data[,2]
  data[,1] = round(data[,1], 3)
  data[,2] = round(data[,2], 3)
  
  # Call the function that will actually create the energy matrix.
  energy_frame = energy_matrix(data, file)
  
  # write the data frame in a file.
  #filename = paste0("coord_lists/", name_prefix, "_s", width, "_coord_list.txt")
  filename = paste0("coord_lists/", name_prefix, "_coord_list.txt")
  write.table(energy_frame, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}
