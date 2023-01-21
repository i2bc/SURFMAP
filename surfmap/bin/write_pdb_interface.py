import argparse
from pathlib import Path

from surfmap.lib.docker import DockerCLI
from surfmap.lib.interface import pdb_from_interface


def get_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-pdb",
        default="str",
        required=True,
        help="Path to a PDB filename"
    )

    parser.add_argument(
        "-res",
        default="str",
        required=True,
        help="Path to a text file informing on interface residues."
    )

    parser.add_argument(
        "-out",
        type=str,
        required=False,
        default=".",
        help = "Output directory"
    )

    parser.add_argument(
        "-suffix",
        type=str,
        required=False,
        default="bs.pdb",
        help = "Suffix added to the written PDB filename. This suffix is added to the original PDB filename stem (e.g. 1g3n_suffix for '-pdb 1g3n.pdb')."
    )

    parser.add_argument(
        "--docker",
        action="store_true",
        help="If chosen, SURFMAP will be run on a docker container (requires docker installed)."
    )

    return parser.parse_args()


def run_from_local(args: argparse.Namespace):
    pdb_filename = args.pdb
    res_filename = args.res
    outdir = Path(args.out)
    outdir.mkdir(exist_ok=True, parents=True)    
    out_filename = outdir / f"{Path(args.pdb).stem}_{args.suffix}"

    pdb_from_interface(pdb_filename=pdb_filename, res_filename=res_filename, out_filename=out_filename)


def run_from_container(args: argparse.Namespace):
    cli = DockerCLI(
        input_args=["-pdb", "-res"],
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