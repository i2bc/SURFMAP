#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, shutil, subprocess


def show_copyrights():
    msg = """
SURFMAP:    Projection of protein surface features on 2D map
Authors:    Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
            Simon Gosset, Anne Lopes
Version:    1.0
Copyright (c) 2021, H. Schweke


SURFMAP relies on the use of MSMS (Copyright (c) M. F. Sanner, 1994).
MSMS is free for academic use. For commercial use please contact
M. F. Sanner at sanner@scripps.edu.

If you use SURFMAP for your research, please cite the following
papers:
    - Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
      Simon Gosset, Anne Lopes (2021) SURFMAP: a software for mapping
      in two dimensions protein surface features. bioRxiv
      2021.10.15.464543; doi: https://doi.org/10.1101/2021.10.15.464543

    - Sanner, M.F., Spehner, J.-C., and Olson, A.J. (1996) Reduced 
      surface: an efficient way to compute molecular surfaces. 
      Biopolymers, Vol. 38., (3), 305-320.

SURFMAP can also optionnaly map electrostatic potential through the use
of the APBS software. So if you use SURFMAP with the option 
'-tomap electrostatics' and find it useful to your research please cite
one or more of the papers listed in the following webpage:
https://apbs.readthedocs.io/en/latest/supporting.html
---------------------------------------------------------------------------

    """
    print(msg)
       
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb",required = True, help = "Input pdb file (path + file name)")
    parser.add_argument("-tomap", type = str, required = True, choices = set(("all", "stickiness", "kyte_doolittle", "wimley_white", "electrostatics", "circular_variance", "circular_variance_atom", "bfactor", "binding_sites")), help = "Choice of the scale. Argument must be one of the following: stickiness; kyte_doolittle; wimley_white; electrostatics; circular_variance; bfactor; binding_sites; all")
    parser.add_argument("-coords", help = argparse.SUPPRESS)
    parser.add_argument("-res", help = "File containing a list of residues to map on the projection. Format must be the following: col 1 = chain id; col 2 = res number; col 3 = res type")
    parser.add_argument("-rad", required = False, help = "Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the smoother the surface (default: 3.0)")
    parser.add_argument("-d", required = False, help = "Output directory where all files will be written (default: './output_SURFMAP_$pdb_$tomap' where $pdb and $tomap are the inputs given to -pdb and -tomap arguments, respectiveley)")
    parser.add_argument("-s", required = False, help = "Size of a grid cell. Necessary that 180%%cellsize == 0 (default: 5.0)")
    parser.add_argument("--nosmooth", action = "store_true", help = "If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)")
    parser.add_argument("--png", action = "store_true", help = "If chosen, a map in png format is computed (default: only pdf format is generated)")
    parser.add_argument("--keep", action = "store_true", help = "If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)")
    args = parser.parse_args()

    pdbarg = args.pdb
    pdb_id = os.path.basename(pdbarg).split(".pdb")[:-1][0]
    pdbname = os.path.basename(pdbarg)
    resfile = args.res
    ppttomap = args.tomap
    
    if not os.path.isfile(pdbarg):
        print("pdb file not found. It seems that the input pdb file does not exist.\nThis could be due to a mistake in the path to the file, for example.\nExiting now.")
        exit()

    if args.res:
        if not os.path.isfile(args.res):
            print("The residue file could not be found (arg -res). It seems that this file does not exist.\nThis could be due to a mistake in the path to the file, for example.\nExiting now.")
            exit()
    
    if args.s:
        try:
            cellsize = int(args.s)
        except ValueError:
            print("cellsize (arg -s) is not an integer.\nYou need to provide an integer if you want to specify the grid resolution.\nExiting now.")
            exit()
        if 180%cellsize != 0:
            print("cellsize (arg -s) is not a multiplier of 180.\nThe cell size needs to be a multiplier of 180 in order to compute the grid. \nExiting now.")
            exit()
        elif int(cellsize) > 20:
            print("the cellsize is too large (arg -s).\nThe maximum cell size authorized is 20, which represents a 9*18 squares grid. \nExiting now.")
            exit()
    else:
        cellsize = 5

    try:
        coordstomap = args.coords
    except:
        coordstomap = None
    
    if args.d:
        outdir = args.d
    else:
        if args.tomap == "all":
            outdir = "output_SURFMAP_"+pdb_id+"_all_properties"
        else:
            outdir = "output_SURFMAP_"+pdb_id+"_"+ppttomap
    
    try:
        if os.path.isdir(outdir) == False:
            os.makedirs(outdir)
    except:
        print("Could not create the output directory. This could be due to a mistake in the input directory name (arg -d).\nExiting now")
        exit()

    try:
        shutil.rmtree(inputdir+"coord_lists")
    except:
        pass
    
    try:
        shutil.rmtree(inputdir+"smoothed_matrices")
    except:
        pass

    try:
        shutil.rmtree(inputdir+"matrices")
    except:
        pass
   
    curdir = os.getcwd()
    absdir = os.path.dirname(os.path.abspath(__file__))
    surftool = absdir+"/tools/SurfmapTools.py"
    shelltool = absdir+"/scripts/compute_shell.sh"
    coordtool = absdir+"/scripts/computeCoordList.R"
    mattool = absdir+"/scripts/computeMatrices.R"
    maptool = absdir+"/scripts/computeMaps.R"


    #============================================================
    # Part 1: Generation of shell around protein surface.
    
    if ppttomap == "electrostatics":
        elecval = "1"
    else:
        elecval = "0"

    if args.rad:
        rad = args.rad
        try:
            float(args.rad)
        except:
            print("The radius given in input (arg -rad) could not be converted to float.\nWhat you provided in input is probably not a number. Please provide a number between 0.5 and 10.\nExiting now.")
            exit()
        if not (float(args.rad) < 10.0 and float(args.rad) > 0.5):
            print("The radius given in input (arg -rad) is either too large or too small. Please provide a number between 0.5 and 10.\nExiting now.")
            exit()
    else:
        rad = "3.0"

    if args.tomap == "all":
        listtomap = ("kyte_doolittle", "stickiness", "wimley_white", "circular_variance")
    else:
        listtomap = []
        listtomap.append(ppttomap)

    # Create shell. 2 options: simple shell or shell with elec potential computed by APBS.
    outputsubp = subprocess.call(["bash", shelltool, "-p", pdbarg, "-e", elecval, "-r", rad])
    shell = curdir + "/shells/" + pdb_id + "_shell.pdb"

    if outputsubp != 0:
        print("the script compute_shell.sh failed. It can be because MSMS could not compute the Connolly surface of the protein.\nYou need to provide a pdb file that MSMS can handle.\nPlease check that there is no problem with your input pdb file.\nExiting now.")
        exit()

    for tomap in listtomap:
        
        print("surface property mapping:",tomap)
        
        BS = False
        
        if tomap == "binding_sites":
            BS = True # keep info of binding site mapping
            tomap = "bfactor"
            scale_opt = "--discrete"
        elif tomap == "circular_variance_atom":
            scale_opt = "--circular_variance"
        else:
            scale_opt = "--" + tomap

        #=============================================================
        # Part 2: convert cartesian coordinates of each particule into spherical coordinates and associates the value of interest (electrostatics, hydrophobicity, stickiness...)

        if not os.path.isfile(shell):
            print("Shell not found. This is probably because MSMS did not manage to compute the surface of the protein.\nExiting now.")
            exit()

        if args.res:
            cmd_surfmap = ["python3", surftool, "-pdb", pdbarg, "-shell", shell, "-tomap", tomap, "-res", resfile, "-d", outdir]
            reslist = os.path.splitext(resfile)[0] + "_sph_coords.out"
        else:
            cmd_surfmap = ["python3", surftool, "-pdb", pdbarg, "-shell", shell, "-tomap", tomap, "-d", outdir]

        outsubres = subprocess.call(cmd_surfmap)
        
        if outsubres != 0:
            print("\nMapping of the residues given in input failed. This can be due to a malformed residue file or a mistake in the reisdue numbering, type or chain.\n\nFormat should be the following:\nchain residuenumber residuetype\nexample: A\t5 LEU\n\nExiting now.")
            exit()


        #=============================================================
        # Part 3: computing phi theta list

        mapfile = "%s_%s_partlist.out"%(pdbname.split(".")[0], tomap)
        cmdlist = ["Rscript", coordtool, "-f", outdir+"/"+mapfile, "-s", str(cellsize)]
        subprocess.call(cmdlist)

        #=============================================================
        # Part 4: computing matrix
        
        coordfile = "%s_%s_coord_list.txt"%(pdbname.split(".")[0], tomap)
        if args.nosmooth:
            if BS:
                cmdmat = ["Rscript", mattool, "-i", outdir+"/coord_lists/"+coordfile, "-s", str(cellsize), "--discrete"]
            else:
                cmdmat = ["Rscript", mattool, "-i", outdir+"/coord_lists/"+coordfile, "-s", str(cellsize), "--nosmooth"]
        else:
            if BS:
                cmdmat = ["Rscript", mattool, "-i", outdir+"/coord_lists/"+coordfile, "-s", str(cellsize), "--discrete"]
            else:
                cmdmat = ["Rscript", mattool, "-i", outdir+"/coord_lists/"+coordfile, "-s", str(cellsize)]
        subprocess.call(cmdmat)
        

        #=============================================================
        # Part 5: computing map

        matfile = "%s_%s_smoothed_matrix.txt"%(pdbname.split(".")[0], tomap)
        matdir = outdir + "/smoothed_matrices/"
        matf = matdir + matfile
        
        matfile2 = "%s_%s_matrix.txt"%(pdbname.split(".")[0], tomap)
        matdir2 = outdir + "/matrices/"
        matf2 = matdir2 + matfile2
        
        if args.png:
            if args.coords:
                if args.res:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "--png", "-c", coordstomap, "-l", reslist, "-s", str(cellsize), "-p", pdb_id]
                else:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "--png", "-c", coordstomap, "-s", str(cellsize), "-p", pdb_id]
            else:
                if args.res:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "--png", "-l", reslist, "-s", str(cellsize), "-p", pdb_id]
                else:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "--png", "-s", str(cellsize), "-p", pdb_id]
        else:
            if args.coords:
                if args.res:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "-c", coordstomap, "-l", reslist, "-s", str(cellsize), "-p", pdb_id]
                else:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "-c", coordstomap, "-s", str(cellsize), "-p", pdb_id]
            else:
                if args.res:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "-l", reslist, "-s", str(cellsize), "-p", pdb_id]
                else:
                    cmdmap = ["Rscript", maptool, "-i", matf, scale_opt, "-s", str(cellsize), "-p", pdb_id]
        subprocess.call(cmdmap)
    
        if not args.keep: # Deleting all intermediate files and directories
            os.remove(outdir + "/" + mapfile)
            os.remove(outdir + "/coord_lists/" + coordfile)
            os.remove(outdir + "/matrices/" + matfile2)
            
            try:
                os.remove(outdir + "/" + pdb_id + "_CV.pdb")
            except OSError:
                pass
            
            if not os.listdir(outdir + "/coord_lists/"): # removing directory if empty
                try:
                    shutil.rmtree(outdir + "/coord_lists/")
                except OSError:
                    pass
            
            if not os.listdir(outdir + "/matrices/"): # removing directory if empty
                try:
                    shutil.rmtree(outdir + "/matrices/")
                except OSError:
                    pass


    # Creating log in output directory. Contains the parameters used to compute the maps.
    outlog = open(outdir+"/log_parameters", "w+")
    outlog.write("parameters used to compute the maps:\n")
    outlog.write("grid resolution: " + str(int(360/cellsize)) + "*" + str(int(180/cellsize)) + "\n")
    outlog.write("MSMS radius: " + rad + "\n")
    outlog.write("property(ies) mapped: " + str(listtomap).strip('[]')+ "\n")
    if args.nosmooth:
        outlog.write("smoothing: off")
    else:
        outlog.write("smoothing: on")

    # Removing shell and electrostatics directories after use. rm residue file in spherical coordinates if exists.
    try:
        shutil.rmtree("shells")
    except OSError:
        pass
    
    try:
        shutil.rmtree("tmp-elec")
    except OSError:
        pass
    
    if args.res:
        try:
            os.remove(reslist)
        except OSError:
            pass

    # Moving the pdb file with CV mapped in bfactor to the output directory
    try:
        shutil.move(os.path.dirname(pdbarg) + "/" + pdb_id + "_CV.pdb", outdir)
    except OSError:
        pass



if __name__ == "__main__":
    show_copyrights()
    main()
