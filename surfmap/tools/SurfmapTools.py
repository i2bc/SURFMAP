#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from pathlib import Path
import sys
from typing import Union

from surfmap.tools import Structure
import numpy


def parseResidueList(reslist, dPDB, outfile):
    # Retrieve CM (cartesian coordinates) of pdb
    (xCM, yCM, zCM) = Structure.centerMassResidueList(dPDB, computeForMulti = True)

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


def run_spherical_coords(shell: Union[str, Path], pdb: Union[str, Path], tomap: str, outdir: Union[str, Path]=".", res: Union[str, Path]=None):
    """Assign atom/residue property value to its closest shell particle and compute the spherical coordinates of this particle.

    Args:
        shell (Union[str, Path]): path to the shell file produced with MSMS
        pdb (Union[str, Path]): path to the input pdb file
        outdir (str): output directory
        tomap (str): property to map
        pdb (Union[str, Path]): path to a residue list to map (optional)
    """
    if not Path(shell).exists():
        print("Could not find the shell file. This is probably because MSMS could not compute it. Generally it is due to a malformed pdb file.\nExiting now")
        exit()

    # extracting information from receptor and particles
    is_bfactor = False if "electrostatics" in tomap else True
    is_charge = True if "electrostatics" in tomap else False

    if "circular_variance" in tomap:
        perres = False if "atom" in tomap else True
        pdb_cv_filename = str(Path(outdir) / f"{Path(pdb).stem}_{tomap}.pdb")
        Structure.compute_CV(pdb, perres=perres, outfilename=pdb_cv_filename)

    input_pdb = pdb_cv_filename if "circular_variance" in tomap else pdb
    dPDB = Structure.parsePDBMultiChains(input_pdb, bfactor=is_bfactor)  # get pdb dict structure
    dshell = Structure.parsePDBParticule(shell, infile=is_charge)  # get shell dict structure

    try:
        Path(pdb_cv_filename).unlink(missing_ok=True)
    except:
        pass
    
    # redefine tomap if necessary
    if tomap in ["interface", "circular_variance", "circular_variance_atom"]:
        property = "bfactor"
    else:
        property = tomap

    # generate specific spherical coords file if some residues are asked to be mapped
    outfile_res_to_map = None
    if res:
        outfile_res_to_map = str(Path(outdir) / f"{Path(res).stem}_sph_coords.out")
        parseResidueList(res, dPDB, outfile=outfile_res_to_map)

    # compute center of mass of receptor
    CMR = Structure.centerMassResidueList(dPDB, all=True, reslist=False, computeForMulti=True, chain=" ")

    # get coords structure to computes distances
    coordlist, idres = Structure.get_coords_idres(dPDB=dPDB)
    
    partlist_outfile = Path(outdir) / f"{Path(pdb).stem}_{tomap}_partlist.out"
    with open(partlist_outfile, "w") as out:
        out.write("phi\ttheta\tvalue\tresnb\trestype\tchain\n")

        # looping over all particules and computing corresponding Value
        for particle in dshell["partlist"]:

            # get coords of the particule
            coord_particle = (dshell[particle]["x"], dshell[particle]["y"], dshell[particle]["z"])  
            _, phi, theta = Structure.coord2spherical(CMR, coord_particle)  # get spherical coordinates from particle cartesian coords 

            # get the closest atom of the shell particle
            _, closestAtom = Structure.getAtomCMRDist(coordlist=coordlist, idres=idres, CMR=coord_particle)     
            chainid, resid, atomtype = closestAtom.split("_")

            # assign the property value of atoms/residues to their closest shell particle
            if property == "electrostatics" :
                scalevalue = dshell[particle]["charge"]
            elif property == "bfactor" :
                scalevalue = dPDB[chainid][resid][atomtype]["bfactor"]
            else:                
                scalevalue = Structure.returnPropensity(aa=dPDB[chainid][resid]["resname"], scale=tomap)            

            out.write("{:8f} {:8f} {:3f} {:8} {:8} {:8}\n".format(phi, theta, scalevalue, dPDB[chainid][resid]["resnum"], dPDB[chainid][resid]["resname"], chainid))

    return outfile_res_to_map, partlist_outfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", required = True, help = "input pdb file (path + file name)")
    parser.add_argument("-shell", required = True, help = "shell file produced with MSMS")
    parser.add_argument("-d", required = True, default=".", help = "output directory")
    parser.add_argument("-res", required = False, default=None, help = "residue list to map (optional)")
    parser.add_argument("-tomap", required = False, type=str, default="electrostatics", help = "value to map (stickiness, kyte-doolittle, wimley-white, bfactor or elec. (Default: elec).")
    parser.add_argument("-surfth", required = False, type=float, default="5.0", help = "output file name (default = map+option in input directory)")
    parser.add_argument("-dist", required = False, help = "distance threshold that separates the spheres from the protein surface (default 5A)")
    parser.add_argument("-onlySurf", required = False, action = "store_true", help = "to use only surface residues for the calculation of the value to map (Default: False).")
    args = parser.parse_args()
    
    shell = args.shell
    pdb = args.pdb
    outdir = Path(args.d)
    tomap = args.tomap
    res=args.res

    run_spherical_coords(shell=shell, pdb=pdb, outdir=outdir, tomap=tomap, res=res)

if __name__ == "__main__":
    main()
