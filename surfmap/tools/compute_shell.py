import argparse
from pathlib import Path
import subprocess
from typing import Union

from surfmap import PATH_MSMS
from surfmap.bin.pdb2xyzr import run as run_pdb2xyzr
from surfmap.bin.multival_csv_to_pdb import run as run_multival_csv2pdb


def get_args():
    parser = argparse.ArgumentParser(description="Script to compute a shell with MSMS from a PDB file. Returns CSV and PDB files with coordinates of the shell surface particles.")

    parser.add_argument(
        "-pdb",
        required=True,
        help="Path to the pdbfile of the protein complex"
    )
    
    parser.add_argument(
        "-rad",
        type=float,
        default=3.0,
        help="Radius added to the standard radius of residues - the higher the radius, the smoother the surface. Defaults to 3.0."    
    )

    parser.add_argument(
        "-out",
        type=str,
        default=Path("."),
        help="Path of the output"
    )

    parser.add_argument(
        "-subdir",
        type=str,
        default=Path("shells"),
        help="Path name of the shell directory content"
    )

    parser.add_argument(
        "--explicit_h",
        action="store_false",
        required=False,
        help="Boolean. use explicit hydrogen radii instead of default united atom radii (default: False)"
    )
    
    return parser.parse_args()


def vert2csv(vertfile: Union[str, Path], outfile: Union[str, Path], skiplines: list=[]):
    """Convert MSMS output file .vert format into .csv format containing only coordinates
    of all the vertices of the surface model created.

    Only the first three columns (x,y,z coordinates) are retrieved and written in a csv format file.
    In MSMS 2.6.1, the first three lines of the .vert file must be skipped.

    The csv file is required as an input for the multivalue executable.

    Args:
        vertfile (Union[str, Path]): Path to the .vert file generated by MSMS
        outfile (Union[str, Path]): Path to the csv output file generated.
        skiplines (list, optional): List of line indexes to skip (index starts at 0). Defaults to [].
    """
    with open(outfile, "w") as _outfile:
        with open(vertfile, "r") as _readfile:
            for i, line in enumerate(_readfile):
                if i in skiplines:
                    continue
                x, y, z = [ _.strip() for _ in line.split()[:3] ]
                _outfile.write(f"{x},{y},{z}\n")



def run(pdb_filename: Union[str, Path], out_dir: Union[str, Path]=".", out_subdir: Union[str, Path]="shells", extra_radius: float=3.0, explicit_hydrogen: bool=False):
    # Create directory containing shells and other files
    outdir = Path(out_dir) / Path(out_subdir)
    outdir.mkdir(exist_ok=True, parents=True)
    outfile_basename = str(outdir / Path(pdb_filename).stem)

    # generate xyzr format input for MSMS
    outfile_xyzr = f"{outfile_basename}.xyzr"
    run_pdb2xyzr(pdb_filename=pdb_filename, out_filename=outfile_xyzr, extra_radius=extra_radius, explicit_hydrogen=explicit_hydrogen)

    # run MSMS to compute Connolly surface
    outfile_vert = f"{outfile_basename}.vert"  # output generated by msms
    # outfile_face = f"{outfile_basename}.face"  # output generated by msms
    msms = f"{str(PATH_MSMS / 'msms')}"

    cmd_msms = [msms, "-if", outfile_xyzr, "-of", f"{outfile_basename}"]
    completed_process = subprocess.run(cmd_msms, capture_output=True)

    # convert .vert file shell coordinates into .csv file
    outfile_csv = f"{outfile_basename}.csv"
    vert2csv(vertfile=outfile_vert, outfile=outfile_csv, skiplines=list(range(3)))

    # convert .csv file into a PDB file
    outfile_pdb = f"{outfile_basename}_shell.pdb"
    run_multival_csv2pdb(ifile=outfile_csv, ofile=outfile_pdb)

    return outfile_csv, outfile_pdb


def main():
    args = get_args()

    return run(
        pdb_filename=args.pdb,
        out_dir=args.out,
        out_subdir=args.subdir,
        extra_radius=args.rad,
        explicit_hydrogen=args.explicit_h
    )

