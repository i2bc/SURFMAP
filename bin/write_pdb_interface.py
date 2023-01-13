from pathlib import Path
from typing import Union

from surfmap.tools import  Structure
 

def parse_res_file(filename: Union[Path, str]) -> dict: 
    """Convert a file listing residues of interest into a dictionary.

    The file is formatted with the following informations (space or tab separated):
    "chain" "resid" "resname"

    Example:
    A  2  1
    A  3  1
    A  6  1

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

            chain, resid, site_label = line.split()
            if chain not in  bs_residues:
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

            try:
                residues_ids, bfactors = zip(*bs_residues[chain])
            except:
                print(f"Error: chain {chain} not found in the PDB structure. Process stopped.")
                exit()

            for atom in structure[chain][resid]["atomlist"]:
                if chain in bs_residues and resid in residues_ids:
                    bfactor = int(bfactors[residues_ids.index(resid)])
                    print(f"Setting bfactor to {bfactor} for atom {atom} of residue {resid} from chain {chain}")
                else:
                    bfactor = 0

                structure[chain][resid][atom]["bfactor"] = bfactor


if __name__ == "__main__":

    try:
        EXAMPLE_PATH = Path(__file__).parent
    except:
        EXAMPLE_PATH = Path(".").parent

    PDB_PATH = EXAMPLE_PATH / "example" / "1g3n_A.pdb"
    RESIDUES_PATH = EXAMPLE_PATH / "example" / "residues_bs.txt"

    
    pdb = Structure.parsePDBMultiChains(infile=PDB_PATH, bfactor=True)
    bs_residues = parse_res_file(RESIDUES_PATH)
    output_filename = EXAMPLE_PATH / (PDB_PATH.stem + "_bs.pdb")


    set_bfactor_from_res(structure=pdb, bs_residues=bs_residues)

    Structure.writePDB(dPDB=pdb, filout=output_filename, bfactor=True)


    # for chain in pdb["chains"]:
    #     for resid in pdb[chain]["reslist"]:
    #         for atom in pdb[chain][resid]["atomlist"]:
    #             print(resid, pdb[chain][resid]["resname"], pdb[chain][resid][atom]["bfactor"])

