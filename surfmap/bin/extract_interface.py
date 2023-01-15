#!/usr/bin/env python

"""
For a given, chain, this script extract the residues that are found 
in the interface with all the other chains of a protein complex.

Interface residues are deduced from ASA difference computations through the use of freesasa.

As freesasa uses PDB file as input, a PDB file must be created for:
- the isolated given chain 
- the other chain(s)

"""
import argparse
from itertools import combinations
from pathlib import Path
import tempfile
from typing import Union

from surfmap.tools import Structure
from surfmap.bin.write_pdb_interface import run as write_pdb_interface
from surfmap.lib.utils import JunkFilePath

import freesasa


def getSurfRes(dASAcp, dASArec):
    """
    Extract surface residues through the comparison of their ASA in the monomer VS dimer
    takes as inputs: dico ASA for cplx (dimer) and rec respectively

    """
        
    #init var
    recsurf = []
    recinter = []

    chainrec = list(dASArec.keys())[0]

    # loop over cplx chains
    for chncp in dASAcp :
        # loop over residues                 
        for rescp in dASAcp[chncp] :
            if chncp == chainrec:

                dASA = dASArec[chncp][str(rescp)].relativeTotal - dASAcp[chncp][str(rescp)].relativeTotal
                    
                if dASA > 0 : # means the ASA in the lig is greater than in the cplx form, so the res is buried in the interface                        
                    recinter.append(rescp)
                    recsurf.append(rescp)
                else: # not in interface
                    if dASArec[chncp][str(rescp)].relativeTotal > 0.25 : # means the res is located on the lig surface
                        recsurf.append(rescp)
            
    return recsurf, recinter


def write_interRes_inBfactor(dPDBrec, chainrec, all_rec_inter,monodir,rootname):
    
    # stores status in bfactor for the studied chain
    nb_interface = len(all_rec_inter)

    for resi in dPDBrec[chainrec]["reslist"] :
        bf = 0
        for i in range(nb_interface):
            interface_id = i+1
            if resi in all_rec_inter[i]:
                bf = interface_id
            for atomi in dPDBrec[chainrec][resi]["atomlist"] :
                dPDBrec[chainrec][resi][atomi]["bfactor"] = bf

    # writes pdb for the studied chain
    Structure.writePDB(dPDBrec, filout="%s/%s_r_inter_bf.pdb"%(monodir, rootname+"_"+chainrec), bfactor=True)
    

def write_interface_residues(interface_map: list, chain: str, outfile: Union[str, Path]) -> None:
    """
    Write a text file listing interface residues of a chain.
    
    The file is formatted with the following informations (space or tab separated):
    "chain" "resid" "resname"

    Example:
    A  2  1
    A  3  1
    A  6  1
    """

    with open(outfile, "w") as _f:
        for i, interface in enumerate(interface_map, start=1):
            for resid in interface:
                _f.write(f"{chain}\t{resid}\t{i}\n")


def generate_dimers(structure: dict, chain_main: str, outdir: str=".") -> list:
    """Generate a set of dimers from a given structure, a chain of interest that will be present in all dimers
    and the other chains to be used.

    For instance, if:
        - structure "xxxx.pdb" has 4 chains: A, B, C and D
        - chain_main is defined as the chain B

    then the script will generate PDB files of the following dimers:
        - AB, BC and BD

    Args:
        structure (dict): Structure of the complex from Structure.parsePDBMultiChains()
        chain_main (str): Main chain that will be present in all dimers
        outdir (str, optional): Path where the temporary output files will be written. Defaults to ".".

    Returns:
        list[str]: List of all the generated dimer filenames
    """
    temp_name = "tmp_dimers_" + next(tempfile._get_candidate_names()) + "_chains-"
    out_basename = str(Path(outdir) / temp_name)

    # get chains in structure that are not the main chain 
    chain_others = [ x for x in structure["chains"] if x != chain_main ]

    # list of all possible dimers that can be formed with the main chain
    dimer_combinations = [ x for x in list(combinations([chain_main] + chain_others , 2)) if chain_main in ''.join(x)]

    # generate dimers
    dimer_filenames = []
    for dimer_chains in dimer_combinations:
        chains_to_rm = [ x for x in structure["chains"] if x not in dimer_chains ]
        dimer_filename = out_basename + ''.join(dimer_chains)

        dimer_filenames.append(dimer_filename)
        Structure.writePDB(dPDB=structure, filout=dimer_filename, bfactor=False, chains_to_rm=chains_to_rm)

    return dimer_filenames


def getResiduesArea(pdbfile: str) -> dict:
    """Get SASA for all residues including relative areas if available for the classifier used.

    Returns 
    Relative areas are normalized to 1, but can be > 1 for residues in unusual conformations or at the ends of chains.

    Args:
        pdbfile (str): Absolut or relative path of a PDB filename

    Returns:
        dict:   dictionary of results where first dimension is chain label and the second dimension residue number:
                >>>result["A"]["5"] 
                -> gives the ResidueArea of residue number 5 in chain A.

    """
    structure = freesasa.Structure(pdbfile)
    sasa = freesasa.calc(structure)
    
    return sasa.residueAreas()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", required=True, help="path to the pdbfile of the protein complex")
    parser.add_argument("-chain", required=True, help="the chain for which this function extract the residues that are found in the interface with all the other chains of a protein complex")
    parser.add_argument("-out", help="path of the output", type=str, default=Path("."))
    parser.add_argument("--write", help="add this argument to create a pdb with the bfactor corresponding to the interfaces", action="store_true")
    args = parser.parse_args()


    junk = JunkFilePath()

    FILES_TO_REMOVE = []


    cplxfile: Union[str, Path] = args.pdb
    chain_to_map: str = args.chain
    outpath: Path = Path(args.out)
    outpath.mkdir(parents=True, exist_ok=True)
    rootname = outpath / Path(cplxfile).stem

    # filenames of pdb files to write for freesasa
    pdb_filename_template = "{}_chain-{}.pdb"
    pdb_chain_filename = pdb_filename_template.format(str(rootname), chain_to_map)

    # get the full pdb structure
    dCPLX = Structure.parsePDBMultiChains(infile=cplxfile)

    # exit program if given chain is absent from pdb file
    if chain_to_map not in dCPLX["chains"]:
        print(f"Error: PDB chain {chain_to_map} not found in the pdb file {cplxfile}")
        exit()

    # list chains different from chain_to_map
    chains_not_rec = [ x for x in dCPLX["chains"] if x != chain_to_map ]

    # write the pdb file for the given chain
    Structure.writePDB(dPDB=dCPLX, filout=pdb_chain_filename, bfactor=False, chains_to_rm=chains_not_rec)
    junk.add(pdb_chain_filename)

    # generate all combinations of possible dimers
    dimer_filenames = generate_dimers(structure=dCPLX, chain_main=chain_to_map, outdir=outpath)
    junk.add(dimer_filenames)

    # compute SASA for the studied chain
    residues_area_chain = getResiduesArea(pdbfile=pdb_chain_filename)

    # compute SASA for each dimers
    residues_area_dimers = []
    for dimer_file in dimer_filenames:
        residues_area_dimers.append(getResiduesArea(pdbfile=dimer_file))

    # compute interface residues of the given chain relative to the other chain in each dimers
    interfaces_chain = []
    for residue_area_dimer in residues_area_dimers:
        _, interface_chain_dimer = getSurfRes(dASAcp=residue_area_dimer, dASArec=residues_area_chain)

        if interface_chain_dimer:
            interfaces_chain.append(interface_chain_dimer)

    
    # write file listing interface residues of the given chain
    interface_filename = f"{str(rootname)}_chain-{chain_to_map}_interface.txt"
    write_interface_residues(interface_map=interfaces_chain, chain=chain_to_map, outfile=interface_filename)


    if args.write:  # write a new PDB file with interface residues information on bfactor column
        out_filename = outpath / f"{Path(pdb_chain_filename).stem}_bs.pdb"
        write_pdb_interface(pdb_filename=pdb_chain_filename, res_filename=interface_filename, out_filename=out_filename)

    # remove intermediray files
    junk.empty()

    
if __name__ == "__main__":
    main()
