
import argparse
from pathlib import Path
import re
from typing import Dict, Tuple, Union

from surfmap import PATH_MSMS
PATH_TO_ATOMTYPENB =  PATH_MSMS / "atmtypenumbers"


def get_radius(atmtype_filename: str) -> Tuple[Dict]:
    explicit_radius = {}
    united_radius = {}
    residue_pattern = {}
    atom_pattern = {}
    atom_number = {}

    with open(atmtype_filename, "r") as atmtype_file:
        for i, line in enumerate(atmtype_file):
            if not line.strip() or line.startswith("#"):
                continue

            if line.split()[0] == "radius":
                atom_key = line.split()[1]
                explicit_radius[atom_key] = line.split()[3]

                if len(line.split()) > 4 and not line.split()[4].startswith("#"):
                    united_radius[atom_key] = line.split()[4]
                else:
                    united_radius[atom_key] = explicit_radius[atom_key]
                continue
            
            residue_pattern[i] = f"^{line.split()[0].replace('*', '.*')}$"
            atom_pattern[i] = f"^{line.split()[1].replace('_', '')}$"
            atom_number[i] = line.split()[2]

            if atom_number[i] not in explicit_radius:
                print(f"Warning: entry {line.strip()} in {atmtype_filename} has no corresponding radius value...")
                print(f"Setting missing radius values to 0.01 for this entry.")
                explicit_radius[atom_number[i]] = 0.01
                united_radius[atom_number[i]] = 0.01

    return explicit_radius, united_radius, residue_pattern, atom_pattern, atom_number


def find_index_pattern(residue_name: str, atom_name: str, residue_pattern: Dict, atom_pattern: Dict) -> Union[int, None]:
    index = None
    for index_pattern, pattern in residue_pattern.items():
        if re.search(pattern, residue_name) and re.search(atom_pattern[index_pattern], atom_name):
            index = index_pattern
            break
    
    return index


def run(pdb_filename: Union[str, Path], out_filename: Union[str, Path], extra_radius: float=3.0, explicit_hydrogen: bool=False):
    # get hard-coded radius information for atoms/residues
    explicit_radius, united_radius, residue_pattern, atom_pattern, atom_number = get_radius(atmtype_filename=PATH_TO_ATOMTYPENB)

    # create output directory if not exist
    Path(out_filename).parent.mkdir(exist_ok=True, parents=True)

    with open(out_filename, "w") as out_file:
        with open(pdb_filename, "r") as _f:
            for line in _f:
                if not line.strip() or line.startswith("#") or line.strip().split()[0].upper() not in ["ATOM", "HETATM"]:
                    continue

                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()

                index_pattern = find_index_pattern(residue_name=residue_name, atom_name=atom_name, residue_pattern=residue_pattern, atom_pattern=atom_pattern)

                if not index_pattern:
                    print(f"Warning: atom pattern was not found for {atom_name} {residue_name}...")
                    print(f"Setting radius to 0.01 for {atom_name} of residue {residue_name}.")
                    out_file.write(f"{x}\t{y}\t{z}\t{0.01}\n")
                    continue

                if not explicit_hydrogen:
                    out_file.write(f"{x}\t{y}\t{z}\t{float(united_radius[atom_number[index_pattern]]) + extra_radius}\n")
                else:
                    out_file.write(f"{x}\t{y}\t{z}\t{float(explicit_radius[atom_number[index_pattern]]) + extra_radius}\n")


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pdb",
        type=str,
        required=True,
        help="Path of a PDB filename"
    )

    parser.add_argument(
        "-rad",
        required=False,
        type=float,
        default=3.0,
        help="Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the smoother the surface (default: 3.0)"
    )

    parser.add_argument(
        "-out",
        required=False,
        type=str,
        default=None,
        help = "Path of the output filename (default: './{pdb_filename}.xyzr')"
    )

    parser.add_argument(
        "--explicit_h",
        action="store_false",
        required=False,
        help="Boolean. use explicit hydrogen radii instead of default united atom radii (default: False)"
    )
    
    return parser.parse_args()


def main():
    args = get_args()

    out_filename = Path(args.out) if args.out else f"{Path(args.pdb).stem}.xyzr"

    run(pdb_filename=args.pdb, out_filename=out_filename, extra_radius=args.rad, explicit_hydrogen=args.explicit_h)
