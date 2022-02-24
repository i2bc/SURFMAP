rm(list=ls())

###################  LIBRARIES   ##########################

library(optparse)
library(mapproj)
options(warn=-1)

###########################################################



################## FUNCTIONS ##########################

# Smooth the cell given in argument. 
# (i.e. average of the values inside a cell, then average with adjacent cells)
smoother_adj <- function(row, datamat) {

  abs = as.double(row[1])/width
  ord = as.double(row[2])/height
  if ((ord != 1) & (ord != stepord) & (abs != 1) & (abs != stepabs)) {
    # Selecting adjacent cells in the matrix (rk: does not correspond to adjacent cells in the dataframe)
    adj_cells = c(datamat[(abs-2)*stepord+ord-1,3],datamat[(abs-2)*stepord+ord,3],
                  datamat[(abs-2)*stepord+ord+1,3],datamat[(abs-1)*stepord+ord-1,3],
                  datamat[(abs-1)*stepord+ord+1,3],datamat[abs*stepord+ord-1,3],
                  datamat[abs*stepord+ord,3],datamat[abs*stepord+ord+1,3])
    if (Inf %in% adj_cells) { # test if cell is in the border of the projection.
      # If the cell is in the border of the projection it does not have 8 neighbors. We consider the cell as outside value.
      smoothed_value = Inf
    } else {
      # Actual smoothing of the cell (ignoring NA values. Rk: for energy maps we considered NA as 0 values).
      seq_case = c(row[3], adj_cells)
      smoothed_value = mean(seq_case, na.rm=T)
    }  
  } else {
    smoothed_value = Inf
  }
}


# Smooth the cell given in argument. version non smoothed (ie cell value is not averaged with neighbouring cells)
smoother_raw <- function(row, datamat) {

  abs = as.double(row[1])/width
  ord = as.double(row[2])/height
  if ((ord != 1) & (ord != stepord) & (abs != 1) & (abs != stepabs)) {
    # Selecting adjacent cells in the matrix (rk: does not correspond to adjacent cells in the dataframe)
    adj_cells = c(datamat[(abs-2)*stepord+ord-1,3],datamat[(abs-2)*stepord+ord,3],
                  datamat[(abs-2)*stepord+ord+1,3],datamat[(abs-1)*stepord+ord-1,3],
                  datamat[(abs-1)*stepord+ord+1,3],datamat[abs*stepord+ord-1,3],
                  datamat[abs*stepord+ord,3],datamat[abs*stepord+ord+1,3])
    if (Inf %in% adj_cells) { # test if cell is in the border of the projection.
      # If the cell is in the border of the projection it does not have 8 neighbors. We consider the cell as outside value.
      smoothed_value = Inf
    } else if (is.na(row[3]) == TRUE) { # The only case where the value of the cell is modified. We assign the most seen value in the neighbors.
      # Actual smoothing of the cell (ignoring NA values. Rk: for energy maps we considered NA as 0 values).
      seq_case = c(row[3], adj_cells)
      smoothed_value = mean(seq_case, na.rm=T)
      # First remove the NA from the vect of value
      x=adj_cells[!is.na(adj_cells)]
      mostseen = as.numeric(names(sort(summary(as.factor(x)), decreasing=T)[1]))
      if (length(mostseen) == 0) {
        smoothed_value = 0
      } else {
        smoothed_value = mostseen
      }
    } else {
      smoothed_value = row[3]
    }
  } else {
    smoothed_value = Inf
  }
}


# Function to find if cell is inside or outside the projection. If cell is outside the projection a value of Inf is assigned to it.
findproj <- function(row) {
  xcoord = as.double(row[1])
  ycoord = as.double(row[2])

  # Selecting the "best" couple of coordinates of the cell (ie the closest to the coordinates (0,0)).
  if (xcoord < 180) {
    if (ycoord < 90) {
      xcoord = xcoord-180
      ycoord = ycoord-90
    } else {
      xcoord = xcoord-180
      ycoord = ycoord-90-height
    }
  } else {
    if (ycoord < 90) {
      xcoord = xcoord-180-width
      ycoord = ycoord-90
    } else {
      xcoord = xcoord-180-width
      ycoord = ycoord-90-height
    }
  }
  ord = 90 - ycoord
  
  # Back to radians
  ord = ycoord * pi/180
  abs = xcoord * pi/180
  
  # Reverse projection
  theta = ord
  phi = abs/cos(theta)
    
  # If phi is not included in (-pi, pi) that means the cell is outside the projection. 
  # (phi E (-pi,pi) and theta E (-pi/2, pi/2) is sufficient to describe all possible coordinate on a sphere). Here theta always belongs to (-pi/2,pi/2) so we only test phi.
  if (findInterval(phi, c(-pi, pi)) != 1){
    row[3] = Inf
  }

  return(as.double(row[3]))
}


findborders <- function(Data) {
  
  # First create blank matrix -> creates dummy projection -> find borders
  # step = 2
  step = 0.5
  
  blankmat <- matrix(nrow = 360/step*180/step, ncol = 3)
  blankmat[,1] <- rep(seq(-180+step,180,by=step),each = 180/step) 
  blankmat[,2] <- rep(seq(-90+step,90,by=step),360/step) 
  blankmat[,3] = 1
  colnames(blankmat) = c("absc", "ord", "val")
  
  projlist <- c("aitoff", "albers", "azequalarea", "azequidist", "bicentric",
                "bonne", "conic", "cylequalarea", "cylindrical", "eisenlohr", "elliptic",
                "fisheye", "gall", "gilbert", "guyou", "harrison", "hex", "homing",
                "lagrange", "lambert", "laue", "lune", "mercator", "mollweide", "newyorker",
                "orthographic", "perspective", "polyconic", "rectangular", "simpleconic",
                "sinusoidal", "tetra", "trapezoidal")
  pf <- c(0, 2, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 2, 0, 2,
          0, 0, 1, 0, 1, 0, 1, 2, 0, 0, 2)
  
  # dummy projection
  if (pf[projlist == projection] == 0){
    blankproj = mapproject(blankmat[,1], blankmat[,2], projection=projection)
  } else if (pf[projlist == projection] == 1) {
    blankproj = mapproject(blankmat[,1], blankmat[,2], projection=projection, 0)
  } else {
    blankproj = mapproject(blankmat[,1], blankmat[,2], projection=projection, c(0,0))
  }
  
  blankproj = data.frame(absc = blankproj$x*180/pi,
                         ord = blankproj$y*180/pi,
                         val = blankmat[,3])

  # If scaling is not good
  proportion1 = summary(blankproj$absc)[6] / summary(blankmat[,1])[6]
  proportion2 = summary(blankproj$ord)[6] / summary(blankmat[,2])[6]
  
  blankproj[,1] = ceiling(blankproj[,1] / (proportion1*width))*width+180
  blankproj[,2] = ceiling(blankproj[,2] / (proportion2*height))*height+90
  
  summary(blankproj)
  plot(blankproj$absc, blankproj$ord, col = "red", pch = 19, cex = 0.1)
  
  inside_vals = paste(Data[,1], Data[,2]) %in% paste(blankproj[,1], blankproj[,2])
  Data$value[!inside_vals] = Inf
  
  return(Data$value)
}



# Create an energy matrix.
comp_val_matrix <- function(Data) {
  Data$absc = as.integer(Data$absc)
  Data$ord = as.integer(Data$ord)
  
  blankmat <- matrix(nrow = 360/width*180/height, ncol = 2)
  blankmat[,1] <- rep(seq(0+width,360,by=width),each = 180/width) 
  blankmat[,2] <- rep(seq(0+height,180,by=height),360/width) 
  colnames(blankmat) = c("absc", "ord")

  # Creating matrix with averaged values in each pixels
  # First splitting input dataframe into list of dataframe per abscissa range.
  absdata = split.data.frame(Data, Data[,1])

  # First splitting dataframe of value per abscissa range into list of dataframe per ordinate range.
  # -> list of list of dataframe, each dataframe containing all value of a given cell.
  n = names(absdata)
  datacells = lapply(n, function(x) {absdata[[x]] = split.data.frame(absdata[[x]], absdata[[x]][,2]); return(absdata[[x]])})
  
  # Calculating mean value for each cell.  
  names(datacells) <- n 
  if (opt$discrete) {
  datameancells = sapply(datacells, function(x) {nord = names(x);
                   lapply(nord, function(y) {x[[y]][,3]=max(x[[y]][,3]); return(x[[y]][1,c(1:3)])})})
  } else {
  datameancells = sapply(datacells, function(x) {nord = names(x);
                   lapply(nord, function(y) {x[[y]][,3]=mean(x[[y]][,3]); return(x[[y]][1,c(1:3)])})})
  }
  datameancells = unlist(datameancells, recursive = FALSE)
  datameancells = do.call(rbind, datameancells)

  # Calculating list of residus present in each cell.
  reslistcells = sapply(datacells, function(x) {
                           nord = names(x);
                           lapply(nord, function(y) {
                                      listres = paste(x[[y]][,5], x[[y]][,4], x[[y]][,6], sep = "_");
                                      listres=unique(listres); 
                                      x[[y]][,3]=paste(listres, collapse = ", "); 
                                      return(x[[y]][1,c(1:3)])})})
  datalistcells = unlist(reslistcells, recursive = FALSE)
  datalistcells = do.call(rbind, datalistcells)
  datalistcells$absc = as.integer(datalistcells$absc)
  datalistcells$ord = as.integer(datalistcells$ord)
  colnames(datalistcells)[3] = "residues"

  datacells = merge(as.data.frame(datameancells), as.data.frame(datalistcells), by=c("absc", "ord"))

  # Creating the matrix.
  filledmat = merge(as.data.frame(blankmat), datacells, by=c("absc", "ord"), all = TRUE)
    
  # Finding all pixels outside projection and attributing a value of Inf to differenciate with residues inside projection.
  if (projection == "sinusoidal"){
    filledmat[,3] = apply(filledmat, 1, findproj)
  } else {
    filledmat[,3] = findborders(filledmat)
  }
  
  
  # Smoothe the matrix.
  if (opt$nosmooth | opt$discrete) {
    filledmat$svalue = apply(filledmat[,c("absc", "ord", "value")], 1, smoother_raw, datamat = filledmat)
  } else {
    filledmat$svalue = apply(filledmat[,c("absc", "ord", "value")], 1, smoother_adj, datamat = filledmat)
  }
  
  filledmat$value = round(as.double(filledmat$value),3)
  filledmat$svalue = round(as.double(filledmat$svalue),3)

  return(filledmat)
}

#####################################################



#####################################################
#                      MAIN
#####################################################

#cat("##############################################\n\nexecuting computeMatrices.R\n\n\n")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=".", 
              help="path to coord files or coord file", metavar="character"),
  make_option(c("-s", "--cellsize"), type="character", default=5, 
              help="dimension of a grid cell. max is 180. please chose a multiple of 180. (default = 5)", metavar="character"),
  make_option(c("--nosmooth"), action="store_true", default="FALSE", 
              help="if TRUE map is not smoothed, if FALSE map is smoothed", metavar="boolean"),
  make_option(c("--discrete"), action="store_true", default="FALSE", 
              help="if TRUE discrete value are used (max value per cell instead of mean)", metavar="boolean"),
  make_option(c("-P", "--projection"), type="character", default="sinusoidal", 
              help="type of projection, must be chosen between: sinusoidal, mollweide", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
projection <<- opt$projection

print(projection)

# Testing if input is a file, a directory, or neither.
if (file_test("-f", opt$input)) {
  path = dirname(normalizePath(opt$input))
  setwd(path)
  dir.create("../maps", showWarnings = FALSE)
  files = c(basename(opt$input))
} else if (file_test("-d", opt$input)) {
  setwd(opt$input)
  dir.create("../maps", showWarnings = FALSE)
  files = list.files(pattern = "\\.txt$")
} else {
  cat("input arg -i is not a directory nor a file:\nexiting now\n\n")
  stop()
}

# Size of the grid.
if (!is.na(opt$cellsize)) {
  width <<- as.integer(opt$cellsize)
} else {
  width <<- 5
}
height <<- width
stepabs <<- 360/width
stepord <<- 180/width

#cat("grid size: ", width, "\n")

if (180%%width != 0) {
  cat("\nError: Grid cell is not a multiple of 180:\nexiting now.\n\n")
  stop()
}

dir.create("../matrices", showWarnings = FALSE)
dir.create("../smoothed_matrices", showWarnings = FALSE)

for (file in (1:length(files))) {
  #prot_name = substr(files[file], 1, 6)
  name_prefix = gsub("_coord_list.txt", "", files[file]) 
  Data = read.table(files[file], fill = TRUE, header = TRUE)
  
  # Call the function that will actually create the energy matrix.
  val_frame = comp_val_matrix(Data)

  # write the data frame in a file.
  filename1 = paste("../smoothed_matrices/" ,name_prefix, "_smoothed_matrix.txt", sep = "")
  filename2 = paste("../matrices/", name_prefix, "_matrix.txt", sep = "")
  write.table(val_frame[,c("absc","ord", "svalue", "residues")], file = filename1, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(val_frame[,c("absc","ord", "value", "residues")], file = filename2, sep = "\t", row.names = FALSE, quote = FALSE)
}
