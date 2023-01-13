#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, sys, os
from surfmap.tools import Structure
import numpy


def parseResidueList(reslist, dPDB):
    # Retrieve CM (cartesian coordinates) of pdb
    (xCM, yCM, zCM) = Structure.centerMassResidueList(dPDB, computeForMulti = True)

    # output file with spherical coordinates of each residue.
    outfile = os.path.splitext(reslist)[0] + "_sph_coords.out"

    outf = open(outfile, "w+")

    # Compute CM (cartesian coordinates) of each residue of the list.
    with open(reslist) as resf:
        for line in resf:
            if len(line.split()) < 3:
                print("error: problem in residue list file (malformed)")
                sys.exit(1)
            else:
                chain = line.split()[0]
                resnumber = line.split()[1]
                restype = line.split()[2]
    
            if chain in dPDB["chains"]:
                if resnumber in dPDB[chain]:
                    if restype in dPDB[chain][resnumber]["resname"]:
                        # Get cartesian coordinates of CM of residue
                        (xCMres, yCMres, zCMres) = Structure.computeResidueCM(dPDB, chain, resnumber)
                        # Get spherical coordinates of CM of residue with respect to CM of receptor
                        (rho, phi, theta) = Structure.coord2spherical((xCM, yCM, zCM), (xCMres, yCMres, zCMres))        
                        linetowrite = line.strip("\n") + "\t" + str(round(rho, 3)) + "\t" + str(round(phi, 3)) + "\t" + str(round(theta, 3)) + "\n"
                        outf.write(linetowrite)
                    else:
                        print("Error in residue list: residue not of the good type.\nPlease follow the numerotation of the pdb given in input.")
                        sys.exit(1)
                else:
                    print("Error in residue list: residue does not exist.\nPlease follow the numerotation of the pdb given in input.")
                    sys.exit(1)
            else:
                print("Error in residue list: chain does not exist.\nPlease follow the numerotation of the pdb given in input.")
                sys.exit(1)

    outf.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", required = True, help = "input pdb file (path + file name)")
    parser.add_argument("-shell", required = True, help = "shell file produced with MSMS")
    parser.add_argument("-d", required = True, help = "output directory")
    parser.add_argument("-res", required = False, help = "residue list to map (optional)")
    parser.add_argument("-tomap", required = False, help = "value to map (stickiness, kyte-doolittle, wimley-white, bfactor or elec. (Default: elec).")
    parser.add_argument("-surfth", required = False, help = "output file name (default = map+option in input directory)")
    parser.add_argument("-dist", required = False, help = "distance threshold that separates the spheres from the protein surface (default 5A)")
    parser.add_argument("-onlySurf", required = False, action = "store_true", help = "to use only surface residues for the calculation of the value to map (Default: False).")
    args = parser.parse_args()
    
    shell = args.shell
    pdb = args.pdb
    pdbPath = "/".join(os.path.abspath(pdb).split("/")[0:-1]) #.join("/")
    pdbName = os.path.basename(pdb)
    outdir = args.d+"/"

    if os.path.isfile(shell) == False:
        print("Could not find the shell file. This is probably because MSMS could not compute it. Generally it is due to a malformed pdb file.\nExiting now")
        exit()
    
    try:
        tomap = args.tomap
    except:    
        tomap = "elec"

    try:
        surfthreshold = args.surfth
    except:    
        surfthreshold = 5.00

    if args.onlySurf:
        surf = True
    else:    
        surf = False


    # INIT DATA: extracting information from REC and PARTS
    #------------------------------------------------------

    out = open(outdir + "%s_%s_partlist.out"%(pdbName.split(".")[0],tomap),"w")
    out.write("phi\ttheta\tvalue\tresnb\trestype\tchain\n")

    if tomap == "bfactor":
        dPDB = Structure.parsePDBMultiChains(pdb, bfactor = True)
    elif tomap == "interface":
        tomap = "bfactor"
        dPDB = Structure.parsePDBMultiChains(pdb, bfactor = True)
    elif tomap == "circular_variance":
        pdb = Structure.compute_CV(pdb, perres = True)
        dPDB = Structure.parsePDBMultiChains(pdb, bfactor = True)
        tomap = "bfactor"
    elif tomap == "circular_variance_atom":
        pdb = Structure.compute_CV(pdb, perres = False)
        dPDB = Structure.parsePDBMultiChains(pdb, bfactor = True)
        tomap = "bfactor"
    elif tomap == "electrostatics":
        dPDB = Structure.parsePDBMultiChains(pdb)
    else: 
        dPDB = Structure.parsePDBMultiChains(pdb)


    # parses PDB of the particules
    if tomap == "electrostatics":
        dshell = Structure.parsePDBParticule(shell, infile = True)
    else:
        dshell = Structure.parsePDBParticule(shell)

    # If list of residues to map provided
    if args.res:
        parseResidueList(args.res, dPDB)

    # computes center of mass of REC
    CMR = Structure.centerMassResidueList(dPDB, all = True, reslist = False, computeForMulti = True, chain = " ")
    
    # For speed purpose: creating
    coordlist = list()
    idres = list()
    for chaini in dPDB["chains"] :
        for resi in dPDB[chaini]["reslist"]:
            for atomi in dPDB[chaini][resi]["atomlist"] :
                xres = dPDB[chaini][resi][atomi]["x"]
                yres = dPDB[chaini][resi][atomi]["y"]
                zres = dPDB[chaini][resi][atomi]["z"]
                coordlist.append((xres,yres,zres))
                idres.append((chaini, resi, atomi))
    coordlist = numpy.asarray(coordlist) # Necessary step for function closest_atom()

    # looping over all particules and computing corresponding Value
    #--------------------------------------------------------------
    cptpart = 0
    for parti in dshell["partlist"] : 

        cptpart = cptpart + 1
        # get coords of the particule
        coordparti = (dshell[parti]["x"], dshell[parti]["y"], dshell[parti]["z"])
        # get the phi/theta corresponding angles
        R, phi, theta = Structure.coord2spherical(CMR, coordparti)

        # computing the value to map
        #---------------------------
        distmini, closestAtom = Structure.getAtomCMRDist(idres, coordlist, coordparti)        
        chainid, resid, atomid = closestAtom.split("_")

        if tomap == "electrostatics" :
            scalevalue = dshell[parti]["charge"] # coulombic energy

        elif tomap == "stickiness" :
            scalevalue = Structure.returnPropensity(dPDB[chainid][resid]["resname"], distmini, scale = "stickiness")

        elif tomap == "bfactor" :
            scalevalue = dPDB[chainid][resid][atomid]["bfactor"]
        
        elif tomap == "kyte_doolittle" :
            scalevalue = Structure.returnPropensity(dPDB[chainid][resid]["resname"], distmini, scale = "kyte_doolittle")

        elif tomap == "wimley_white" :
            scalevalue = Structure.returnPropensity(dPDB[chainid][resid]["resname"], distmini, scale = "wimley_white")
     
        out.write("{:8f} {:8f} {:3f} {:8} {:8} {:8}\n".format(phi, theta, scalevalue, dPDB[chainid][resid]["resnum"], dPDB[chainid][resid]["resname"], chainid))
    out.close()




if __name__ == "__main__":
    main()
