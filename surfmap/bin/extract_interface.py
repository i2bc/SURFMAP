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
from pathlib import Path
from typing import Union


from surfmap.lib.interface import extract_interface
from surfmap.lib.docker import DockerCLI


def get_args() -> argparse.Namespace:

    description = """
Program allowing to compute interface residues of (a) chain(s) belonging to a given multimeric protein PDB file.
The default output is a text file listing all residues found at the interface between the given chain(s) and
the other chain(s) in the multimeric protein and a PDB file made of the given chain(s). In this PDB file,
the b-factor column is filled with a number for each different interface found.
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-pdb", required=True, help="Path to the pdbfile of the protein complex")
    parser.add_argument("-chains", required=True, help="Chain(s) for which this function extract the residues that are found in the interface with all the other chains of a protein complex. (e.g. '-chains A' or '-chains AD' ")
    parser.add_argument("-out", help="Path of the output", type=str, default=Path("."))
    parser.add_argument("--docker", action="store_true", help="If chosen, SURFMAP will be run on a docker container (requires docker installed).")
    # parser.add_argument("--write", help="Add this argument to create a pdb with the bfactor corresponding to labelled interfaces", action="store_true")
    
    return parser.parse_args()


def run_from_local(args: argparse.Namespace):
    cplxfile: Union[str, Path] = args.pdb
    chain_to_map: str = args.chains
    outpath: Path = Path(args.out)

    extract_interface(cplxfile=cplxfile, chain_to_map=chain_to_map, outpath=outpath)


def run_from_container(args: argparse.Namespace):
    cli = DockerCLI(
        input_args=["-pdb"],
        output_args="-out",
        out_dirname=args.out,
        args=args
    )
    
    cli.show()
    cli.run()


def main():
    args = get_args()

    if args.docker:
        run_from_container(args)
    else:
        run_from_local(args)


if __name__ == "__main__":
    main()
