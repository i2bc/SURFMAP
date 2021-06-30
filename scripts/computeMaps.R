rm(list=ls())

###################  LIBRARIES   ##########################

library(optparse)
options(warn=-1)

###########################################################


###################  PARAM   ##########################

# if true all energy map written in a unique file. Else one file by map.
one_file = TRUE

################## END PARAM ##########################


################## FUNCTIONS ##########################

# This function is a better version of the function image(). It permits to give color to cell that have NA, NaN and 
# values which are outside the color scale. Here the color white if attributed to all the cell with a value of 100, 
# which correspond to the cells outside of the sinusoidal projection.

image.nan.better <- function(z,  zlim = c(0,-minval), col, na.color='gray', outside.below.color='black', outside.above.color='blue',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  
  image(z=z,  zlim=zlim, col=col, ...) # we finally call image(...)
}


#This function creates a color scale for use with e.g. the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
#Depending on the orientation, x- or y-limits may be defined that
#are different from the z-limits and will reduce the range of
#colors displayed.

image.scale <- function(z, zlim, col = colorScale2, scalename,
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", 
       xlab=scale_main, adj = 0.4, ylab="", line=1.5, cex.lab=0.9, cex.names=0.6, ...) 
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

#####################################################



#####################################################
#                       MAIN
#####################################################

option_list = list(
  make_option(c("-i", "--input"), type="character", default=".", 
              help="path to matrix files or matrix file", metavar="character"),
  make_option(c("-p", "--pdb"), type="character", default=".", 
              help="pdb file name (without .pdb extension)", metavar="character"),
  make_option(c("-c", "--coord"), type="character", default=NA, 
              help="File containing (phi, theta) coordinates to map", metavar="character"),
  make_option(c("-l", "--reslist"), type="character", default=NA, 
              help="File containing coordinates of residues to map", metavar="character"),
  make_option(c("-s", "--cellsize"), type="character", default=5, 
              help="dimension of a grid cell. max is 180. please chose a multiple of 180. (default = 5)", metavar="character"),
  make_option(c("--png"), action="store_true", default = "FALSE",
              help="output file in png format (default = pdf)", metavar="character"),
  make_option(c("--electrostatics"), action="store_true", default = "FALSE",
              help="use electrostatics scale", metavar="character"),
  make_option(c("--stickiness"), action="store_true", default = "FALSE",
              help="use stickiness scale", metavar="character"),
  make_option(c("--kyte_doolittle"), action="store_true", default = "FALSE",
              help="use hydrophobicity scale Kyte-Doolittle", metavar="character"),
  make_option(c("--wimley_white"), action="store_true", default = "FALSE",
              help="use hydrophobicity scale Wimley-White", metavar="character"),
  make_option(c("--circular_variance"), action="store_true", default = "FALSE",
              help="use circular variance scale", metavar="character"),
  make_option(c("--bfactor"), action="store_true", default = "FALSE",
              help="use bfactor scale", metavar="character"),
  make_option(c("--energy"), action="store_true", default = "FALSE",
              help="use energy scale", metavar="character"),
  make_option(c("--whiteblue"), action="store_true", default = "FALSE",
              help="use white to blue scale", metavar="character"),
  make_option(c("--discrete"), action="store_true", default = "FALSE",
              help="use discrete scale (up to 10 colors)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
pdb_id = opt$pdb

if (file_test("-f", opt$input)) {
  path = dirname(normalizePath(opt$input))
  files = c(basename(opt$input))
} else if (file_test("-d", opt$input)) {
  path = opt$input
  files = list.files(path, pattern = "\\matrix.txt$")
} else {
  cat("input arg -i is not a directory nor a file:\nexiting now\n\n")
  stop()
}

minval=0
maxval=0

if (!is.na(opt$coord)) {
  theta_phi_nat = read.table(opt$coord, stringsAsFactors=F)
  col_names = c("rec", "theta", "phi")
  colnames(theta_phi_nat) = col_names
}

if (!is.na(opt$reslist)) {
  res_to_map = read.table(opt$reslist)
  col_names = c("chain", "pos", "type", "rho", "phi", "theta")
  colnames(res_to_map) = col_names
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

if (180%%width != 0) {
  cat("\nError: Grid cell is not a multiple of 180. Exiting now.\n")
  stop()
}

# Moving to working directory
setwd(path)
dir.create("../maps", showWarnings = FALSE)

for (file in (1:length(files))) {
  if (opt$png == TRUE) {
    png(paste("../maps/",gsub("_smoothed_matrix.txt", "", files[file]),"_map.png", 
        sep = ""), res = 300, width = 17.78, height = 17.78, units = "cm")
  } else {
    pdf(paste0("../maps/",gsub("_smoothed_matrix.txt", "", files[file]),"_map.pdf"))
  }
  prot_name = gsub("_smoothed_matrix.txt", "", files[file])
  
  # Retieving coordinates of points to plot (if user has enterd coord file).
  if (!is.na(opt$reslist)) {
    legend = do.call(paste, res_to_map[c(3,2)])
    res_nb = res_to_map[,2]
    res_type = res_to_map[,3]
    # Coordinates of the points to plot.
    res_phi = res_to_map[,5]
    res_theta = res_to_map[,6]
    # Sinusoidal projection of the coordinates.
    res_phi = res_phi-pi
    res_theta = res_theta-pi/2
    res_phi = res_phi*cos(res_theta)
    # Conversion of the coordinates from radians to degrees.
    res_phi = res_phi*180/pi+180
    res_theta = res_theta*180/pi+90
  }
  
  # Retieving coordinates of points to plot (if user has enterd coord file).
  if (!is.na(opt$coord)) {
    filename = tools::file_path_sans_ext(file)
    filename = basename(filename)
    rec = substr(files[file],0,17)
    
    # Coordinates of the points to plot.

    # TO COMMENT IN FINAL VERSION !!!!
    #nat_phi=theta_phi_nat[which(theta_phi_nat$rec==rec),3]
    #nat_theta=theta_phi_nat[which(theta_phi_nat$rec==rec),2]
    nat_phi=theta_phi_nat[which(theta_phi_nat[,1]==rec),3]
    nat_theta=theta_phi_nat[which(theta_phi_nat[,1]==rec),2]
    #cat("nat phi: ", nat_phi, "\n")
    #########################

    # TO UNCOMMENT IN FINAL VERSION !!!!!
    #nat_phi=theta_phi_nat[,4]
    #nat_theta=theta_phi_nat[,3]
    # Sinusoidal projection of the coordinates.

    nat_phi=nat_phi-pi
    nat_theta=nat_theta-pi/2
    nat_phi=nat_phi*cos(nat_theta)
    # Conversion of the coordinates from radians to degrees.
    nat_phi=nat_phi*180/pi+180
    nat_theta=nat_theta*180/pi+90
  }

  # Read the current global file.
  data_matrix = read.table(files[file], fill = TRUE, header = TRUE)
  data_matrix[is.na(data_matrix)] <- 0
  val_matrix = matrix(data_matrix[,3],ncol = stepabs, nrow = stepord, byrow=FALSE)

  # The min value of the matrix correspond to the lowest energy value. Used for the scale
  # in the plot, which correspond to (0,-minE).
  proj=which(val_matrix!=100)

  if (opt$electrostatics == TRUE) { # electrostatics scale
    main_scale = "electrostatic potential distribution scale"
    minval=min(val_matrix)
    maxval=max(val_matrix[proj])
    if (abs(minval) > abs(maxval)) {
      maxval = abs(minval)
    } else {
      minval=-(abs(maxval))
    }
    range=abs(minval-maxval)
    scale_main = paste0("electrostatic","\n","potential")
    main_title = paste0("electrostatic potential map", "\n", pdb_id)
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    scale_at = c(ceiling(minval*100000)/100000, round(maxval-range*5/6,5), round(maxval-range*4/6,5), round(maxval-range*3/6,5), round(maxval-range*2/6,5), round(maxval-range*1/6,5), floor(maxval*100000)/100000)
    colorScale <- colorRampPalette(c("red", "white", "blue"))(1000)

  } else if (opt$stickiness == TRUE) { # stickiness scale
    minval = -1.273
    maxval = 1.273
    range = abs(minval-maxval)
    scale_main = "interaction\npropensity"
    main_title = paste0("stickiness map\n", pdb_id)
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    colorScale <- colorRampPalette(c("blue", "white", "red"))(1000)
    scale_at = c(-1.273,-0.85,-0.43,0,0.43,0.85,1.273)

  } else if (opt$kyte_doolittle == TRUE) { # hydrophobicity scale
    main_scale = "hydrophobicity scale Kyte-Doolittle"
    main_title = paste0("Kyte-Doolittle hydrophobicity map\n", pdb_id)
    minval = -4.5
    maxval = 4.5
    range=abs(minval-maxval)
    scale_main = paste0("hydrophobicity","\n","Kyte-Doolittle")
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    colorScale <- colorRampPalette(c("blue", "white", "red"))(1000)
    scale_at = c(-4.5,-3,-1.5,0,1.5,3,4.5)
    val_matrix = round(val_matrix,digits=2)

  } else if (opt$wimley_white == TRUE) { # hydrophobicity scale
    main_scale = "hydrophobicity scale Wimley-White"
    main_title = paste0("Wimley-White hydrophobicity map\n", pdb_id)
    minval = 2.23
    maxval = 6.1
    range=abs(minval-maxval)
    scale_main = paste0("hydrophobicity","\n","Wimley-White")
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    colorScale <- colorRampPalette(c("blue", "white", "red"))(1000)
    scale_at = c(2.23,2.85,3.5,4.15,4.8,5.45,6.1)

  } else if (opt$circular_variance == TRUE) { # circular variance scale
    main_scale = "CV"
    main_title = paste0("circular variance map\n", pdb_id)
    minval = 0
    maxval = 1
    range=abs(minval-maxval)
    scale_main = paste0("circular","\n","variance")
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    colorScale <- colorRampPalette(c("red", "white", "blue"))(1000)
    scale_at = c(0,1/6,2/6,3/6,4/6,5/6,1)

  } else if (opt$bfactor == TRUE) { # b-factor scale
    main_scale = "b-factor"
    main_title = paste0("b-factor map\n", pdb_id)
    
    minval = min(val_matrix[proj])
    maxval = max(val_matrix[proj])
    range = abs(minval-maxval)

    scale_main = "b-factor"
    colors = c(seq(minval,minval+range*1/3,length=334),seq(minval+range*1/3,minval+range*2/3,length=333),seq(minval+range*2/3,maxval,length=334))
    colorScale <- colorRampPalette(c("blue", "white", "red"))(1000)
    scale_at = c(ceiling(minval*100000)/100000, round(maxval-range*5/6,5), round(maxval-range*4/6,5), round(maxval-range*3/6,5), round(maxval-range*2/6,5), round(maxval-range*1/6,5), floor(maxval*100000)/100000)
  
  } else if (!is.na(opt$discrete)) { # b-factor scale
    main_scale = "b-factor"
    main_title = paste0("b-factor map\n", prot_name)
    
    minval = min(val_matrix[proj])
    maxval = max(val_matrix[proj])
    range = abs(minval-maxval)

    scale_main = "b-factor"
    colors = seq(0,maxval+1,1)
    colorScaletot = c("white", "firebrick2", "cornflowerblue", "darkgoldenrod1", "mediumorchid4", "green3", "skyblue1", "tan3", "violet", "darkgreen", "chocolate1", "gray30", "salmon", "palevioletred3", "royalblue4")
    colorScale <- colorScaletot[seq(1,maxval+1,1)]
    scale_at = seq(0.5, maxval+0.5,1)


  } else if (opt$whiteblue == TRUE) { # white to blue scale
    main_scale = "bfactor"
    main_title = paste0("bfactor map\n", pdb_id)
    minval = min(val_matrix[proj])
    maxval = max(val_matrix[proj])
    range = abs(minval-maxval)
    scale_main = paste0("b-factor")
    colors = c(seq(minval,minval+range*1/2,length=500),seq(minval+range*1/2,maxval,length=501))
    colorScale <- colorRampPalette(c("white", "blue"))(1000)
    scale_at = c(0,1/5,2/5,3/5,4/5,1)

  } else { # energy scale
    main_scale = "ATTRACT energy scale"
    main_title = paste0("ATTRACT docking energy landscape\n", prot_name, " - ", lig_name)
    minval = min(val_matrix[proj])
    maxval = max(val_matrix[proj])

    scale_main = expression(paste("   kcal.",mol^-1))
    range=abs(minval-maxval)
    
    colors = c(seq(minval,minval+range*1/4,length=250),seq(minval+range*1/4,minval+range*2/4,length=250),seq(minval+range*2/4,minval+range*3/4,length=250),seq(minval+range*3/4,maxval,length=251))
    colorScale <- colorRampPalette(c("red","yellow","green3","blue"))(1000)
    scale_at = c(ceiling(minval*100000)/100000, round(maxval-range*5/6,5), round(maxval-range*4/6,5), round(maxval-range*3/6,5), round(maxval-range*2/6,5), round(maxval-range*1/6,5), floor(maxval*100000)/100000)

  }

  # Creation of the map.
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1,1), heights=c(4,1))
  par(mar=c(14.4,5,9.3,1.2))
  
  # compute values for scale axis
  range=maxval-minval
  if (opt$energy) { 
      axis_scale = c(ceiling(minval*100000)/100000, round(maxval-range*5/6,5), round(maxval-range*4/6,5), round(maxval-range*3/6,5), round(maxval-range*2/6,5), round(maxval-range*1/6,5), floor(maxval*100000)/100000)
  } else if (opt$whiteblue) { 
      axis_scale = c(ceiling(minval*100000)/100000, round(maxval-range*4/5,5), round(maxval-range*3/5,5), round(maxval-range*2/5,5), round(maxval-range*1/5,5), floor(maxval*100000)/100000)
  } else if (opt$discrete) {
      axis_scale = seq(0,maxval,1)
  } else { 
      axis_scale = c(ceiling(minval*100000)/100000, round(maxval-range*5/6,5), round(maxval-range*4/6,5), round(maxval-range*3/6,5), round(maxval-range*2/6,5), round(maxval-range*1/6,5), floor(maxval*100000)/100000)}

  # Reverse the matrix for correct display.
  val_matrix = apply(val_matrix,2,rev)
  
  # Creation of the image.
  image.nan.better(t(val_matrix),col=colorScale,
                   zlim=c(minval,maxval),
                   outside.below.color='white',
                   outside.above.color='gray90',
                   na.color='white',
                   frame.plot = TRUE,
                   axes = FALSE,
                   xlab = expression(paste(phi, " sin(", theta, ")")), 
                   ylab = expression(paste("90 - ", theta)),
                   cex.lab = 1.5)
  axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(-180,-90,0,90,180), cex.axis=1.2)
  axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(-90,-45,0,45,90), cex.axis=1.2, las=2)
  title(main = main_title, line = 1.5)
  
  # Add native site to plot.
  if (!is.na(opt$coord)) {
    #points(nat_phi/360, nat_theta/180, pch = c(8,13,10,7,3,4), col = "black", cex = 4, lwd = 3)
    points(nat_phi/360, nat_theta/180, pch = 8, col = "black", cex = 3, lwd = 2)
    #text(nat_phi/360, nat_theta/180+22/180, labels = theta_phi_nat[,2], cex = 1.8)
  }

  # Add coordinates of residues of interest
  if (!is.na(opt$reslist)) {
    points(res_phi/360, res_theta/180, pch = 8, col = "black", cex = 2, lwd = 1.5)
    text(res_phi/360, res_theta/180+14/180, labels = legend, cex = 0.7, lwd = 1.5)
  }
  # Add scale to plot. 
  par(mar=c(15.5,1.6,10.5,4.5))
  image.scale(t(val_matrix), col=colorScale, breaks=colors, scalename = scale_main, horiz=FALSE, yaxt="n")#, ylim = c(0, 1))
  axis(4,at=scale_at, las=2, cex.axis=0.8, labels=round(axis_scale,digits = 2))

  dev.off()
}

