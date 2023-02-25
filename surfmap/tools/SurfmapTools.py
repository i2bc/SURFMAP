#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import sys
from typing import Union

from surfmap.tools import Structure
from surfmap.lib.logs import get_logger


logger = get_logger(name=__name__)


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
                        logger.error("Error in residue list: residue not of the good type.\nPlease follow the numerotation of the pdb given in input.")
                        sys.exit(1)
                else:
                    logger.error("Error in residue list: residue does not exist.\nPlease follow the numerotation of the pdb given in input.")
                    sys.exit(1)
            else:
                logger.error("Error in residue list: chain does not exist.\nPlease follow the numerotation of the pdb given in input.")
                sys.exit(1)

    outf.close()


def write_shell_pdb(dshell: dict, outfile: str):
    with open(outfile, "w") as out_shell:
        for i, particle in enumerate(dshell["partlist"], start=1):
            line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                "ATOM",
                i,
                "O",
                "",
                "POS",
                "A",
                1,
                "",
                dshell[particle]["x"],
                dshell[particle]["y"],
                dshell[particle]["z"],
                1.00,
                dshell[particle]["charge"]
            )
            out_shell.write(line)


def run_particles_mapping(shell: Union[str, Path], pdb: Union[str, Path], tomap: str, outdir: Union[str, Path]=".", res: Union[str, Path]=None):
    """Assign atom/residue property value to its closest shell particle and compute the spherical coordinates of this particle.

    Two PDB structures are required:
        - one with coordinates of the shell particles
        - one with coordinates of the initial atoms of a protein

    Each shell particle must be assigned a property value. With the exception of the electrostatics property, each value assigned to a particle is retrieved from the property value
    of its closest atom/residue of the protein. For electrostatics property, values are already assigned to each particle.

    For hydrophobic properties (stickiness, kyte_doolittle and wimley_white), the values are computed on the fly for the closest atom/residue of a given particle.
    For properties like interface, circular_variance or circular_variance_atom, bfactor, the values are retrieved in the bfactor column. So, when circular variance property is asked,
    a new PDB is first generated.

    If a file with a list of residues to map is given, an output file with '_sph_coords.out' as suffix/extension is generated. This file could then be given to the 
    compute_map function to map those residues.

    Args:
        shell (Union[str, Path]): path to the shell file produced with MSMS
        pdb (Union[str, Path]): Path to the input pdb file
        tomap (str): Property to map
        outdir (str): Path to the output directory. Defaults to '.'.
        res (Union[str, Path]): Optional path to a file with list of residues to map (optional)

    Returns:
        outfile_res_to_map (str | None): path to the spherical coords file of the list of residues to be mapped. None if no list of residues exists.
        partlist_outfile (str): Path to the file with spherical coords of every shell particles. The file has the following format: #phi #theta #scalevalue, #resnb, #resname, #chain

    """
    if not Path(shell).exists():
        print("Could not find the shell file. Exiting now\n")
        exit()

    # extracting information from receptor and particles
    is_bfactor = False if "electrostatics" in tomap else True
    is_charge = True if "electrostatics" in tomap else False

    if "circular_variance" in tomap:
        perres = False if "atom" in tomap else True
        pdb_cv_filename = str(Path(outdir) / f"{Path(pdb).stem}_{tomap}.pdb")
        logger.debug(f"Running circular variance computation of {pdb_cv_filename}")
        Structure.compute_CV(pdb, perres=perres, outfilename=pdb_cv_filename)

    # get pdb dict structure of input pdb and shell particles
    input_pdb = pdb_cv_filename if "circular_variance" in tomap else pdb

    logger.debug(f"Get dictionary of the PDB structure from {input_pdb}")
    dPDB = Structure.parsePDBMultiChains(input_pdb, bfactor=is_bfactor)

    logger.debug(f"Get dictionary of the shell structure from {shell}")
    dshell = Structure.parsePDBParticule(shell, infile=is_charge)
    

    # remove pdb of circular_variance if exists
    try:
        Path(pdb_cv_filename).unlink(missing_ok=True)
        logger.debug(f"{pdb_cv_filename} has been removed")
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
        logger.debug(f"Parsing user-given residues to map and generating {outfile_res_to_map}")
        parseResidueList(res, dPDB, outfile=outfile_res_to_map)

    # compute center of mass of receptor
    logger.debug(f"Retrieving the center of mass of the pdb structure")
    CMR = Structure.centerMassResidueList(dPDB, all=True, reslist=False, computeForMulti=True, chain=" ")

    # get coords structure to compute distances
    coordlist, idres = Structure.get_coords_idres(dPDB=dPDB)
    
    partlist_outfile = Path(outdir) / f"{Path(pdb).stem}_{tomap}_partlist.out"
    with open(partlist_outfile, "w") as out:
        out.write(f"{'phi'}\t{'theta'}\t{'value'}\t{'resnb'}\t{'restype'}\t{'chain'}\n")

        # looping over all particules and computing corresponding Value
        logger.debug(f"Looping over all shell particles to assign them the property value of their closest atoms/residues")
        for i, particle in enumerate(dshell["partlist"], start=1):

            # get coords of the particule
            logger.debug(f"Reading coordinates of the particle {i} and converting it into spherical coordinates")
            coord_particle = (dshell[particle]["x"], dshell[particle]["y"], dshell[particle]["z"])
            _, phi, theta = Structure.coord2spherical(CMR, coord_particle)

            # get the closest atom of the shell particle
            logger.debug(f"Retrieving the closest atom of the particle {i}")
            _, closestAtom = Structure.getAtomCMRDist(coordlist=coordlist, idres=idres, CMR=coord_particle)
            chainid, resid, atomtype = closestAtom.split("_")

            # assign the property value of atoms/residues to their closest shell particle or directly from the shell particle
            if property == "electrostatics":
                scalevalue = dshell[particle]["charge"]
                logger.debug(f"Electrostatic value {scalevalue} read from the shell structure has been assigned to the particle {i}")
            elif property == "bfactor":
                    scalevalue = dPDB[chainid][resid][atomtype]["bfactor"]
                    logger.debug(f"{tomap} value {scalevalue} read from the atom {chainid}-{resid}-{atomtype} of {input_pdb} has been assigned to the particle {i}")
            else:
                scalevalue = Structure.returnPropensity(aa=dPDB[chainid][resid]["resname"], scale=tomap)
                logger.debug(f"{tomap} value {scalevalue} computed for the residue {chainid}-{resid} of {input_pdb} has been assigned to the particle {i}")

            # temporary - assign value to shell
            dshell[particle]["charge"] = scalevalue

            out.write("{:8f} {:8f} {:3f} {:8} {:8} {:8}\n".format(phi, theta, scalevalue, dPDB[chainid][resid]["resnum"], dPDB[chainid][resid]["resname"], chainid))

    # temporary - write shell pdb with assigned values
    # shell_pdb = Path(outdir) / f"{Path(pdb).stem}_{tomap}_shell.pdb"
    # write_shell_pdb(dshell=dshell, outfile=shell_pdb)
    
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

    run_particles_mapping(shell=shell, pdb=pdb, outdir=outdir, tomap=tomap, res=res)

if __name__ == "__main__":
    main()
