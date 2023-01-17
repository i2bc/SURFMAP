import argparse
from pathlib import Path
from typing import Union

from surfmap.tools import  Structure
 

def parse_res_file(filename: Union[Path, str]) -> dict: 
    """Convert a file listing residues of interest into a dictionary.

    The file is formatted with the following informations (space or tab separated):
    "chain" "resid" "resname"

    Example:
    #chain   #resid   #resname #label_value
    A   2    PHE    1  
    A   3    GLU    1
    A   6    ALA    1

    The function will return a dictionary as follows:
    {
        "A": [(2, 1), (3, 1), (6, 1)]
    }


    Args:
        filename (Union[Path, str]): Path to the file listing the resiudes of interest.

    Return:
        dict
    """
    bs_residues = {}

    with open(filename, "r") as _file:
        for line in _file:
            if line.startswith("#") or not line:
                continue

            chain, resid, resname, site_label = line.split()
            if chain not in bs_residues:
                bs_residues[chain] = []
            bs_residues[chain].append((resid, site_label))

    return bs_residues                


def set_bfactor_from_res(structure: dict, bs_residues: dict):
    """_summary_

    Args:
        structure (dict): PDB structure - from Structure.parsePDBMultiChains()
        residues (dict): list of residues to set - from parse_res_file()
    """
    for chain in structure["chains"]:
        for resid in structure[chain]["reslist"]:

            if chain in bs_residues:
                residues_ids, bfactors = zip(*bs_residues[chain])

            for atom in structure[chain][resid]["atomlist"]:
                if chain in bs_residues and resid in residues_ids:
                    bfactor = int(bfactors[residues_ids.index(resid)])
                else:
                    bfactor = 0

                structure[chain][resid][atom]["bfactor"] = bfactor


def get_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-p", "--pdb",
        default="str",
        required=True,
        help="Path to a PDB filename"
    )

    parser.add_argument(
        "-r", "--res",
        default="str",
        required=True,
        help="Path to a text file informing on interface residues."
    )

    parser.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        default=".",
        help = "Output directory"
    )

    parser.add_argument(
        "-s", "--suffix",
        type=str,
        required=False,
        default="bs.pdb",
        help = "Suffix added to the written PDB filename. This suffix is added to the original PDB filename stem (e.g. 1g3n_suffix for '-pdb 1g3n.pdb')."
    )

    return parser.parse_args()


def run(pdb_filename: str, res_filename: str, out_filename: str):
    """From a PDB file and a text file informing about interface residues,
    write a PDB file with interface residues information on bfactor.

    This information is 0 for non-interface residues and 1, 2, etc... for residues
    of interface 1, 2, etc..., respectively.

    Args:
        pdb_filename (str): Path to a PDB filename
        res_filename (str): Path to a text file informing on interface residues.
        out_filename (str): Filename of the written PDB file.
    """
    # get pdb structure dict
    pdb = Structure.parsePDBMultiChains(infile=pdb_filename, bfactor=True)

    # parse residues interface file
    interface_residues = parse_res_file(res_filename)

    # edit bfactor in the pdb structure: an integer label value for given interface residues, 0 otherwise
    set_bfactor_from_res(structure=pdb, bs_residues=interface_residues)

    # write the edited struture in pdb format
    Structure.writePDB(dPDB=pdb, filout=out_filename, bfactor=True)


def main():
    args = get_args()

    pdb_filename = args.pdb
    res_filename = args.res

    outdir = Path(args.out)
    outdir.mkdir(exist_ok=True, parents=True)    
    out_filename = outdir / f"{Path(args.pdb).stem}_{args.suffix}"

    run(pdb_filename=pdb_filename, res_filename=res_filename, out_filename=out_filename)


if __name__ == "__main__":
    main()