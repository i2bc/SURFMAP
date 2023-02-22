import argparse
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from typing import Union

from surfmap.bin.multival_csv_to_pdb import run as run_multival_csv2pdb
from surfmap.lib.logs import get_logger


logger = get_logger(name=__name__)


def get_args():
    parser = argparse.ArgumentParser(description="Script to compute electrostatics. Returns PDB files with coordinates of the shell surface particles and electrostatic.")

    parser.add_argument(
        "-pdb",
        required=True,
        type=str,
        help="Path to the pdbfile of the protein complex"
    )
    
    parser.add_argument(
        "-csv",
        required=True,
        type=str,
        help="CSV file with coordinates of the shell particles."    
    )

    parser.add_argument(
        "-pqr",
        required=False,
        type=str,
        default="",
        help="PQR file used to build input file of APBS"    
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
        default=Path("tmp-elec"),
        help="Path name of the electrostatics directory content"
    )
    
    return parser.parse_args()


def edit_inputgen(inputgen_file: Union[str, Path], pqrfile: Union[str, Path]):
    with open(inputgen_file, 'r') as _file :
        data = _file.read()

    # Replace some data
    data = data.replace(Path(pqrfile).name, pqrfile)
    data = data.replace("pot dx pot", f"pot dx {pqrfile}")

    # Write the file out
    with open(inputgen_file, 'w') as _file:
        _file.write(data)


def run(pdb_filename: Union[str, Path], csv_filename: Union[str, Path], pqr_filename: Union[str, Path]="", out_dir: Union[str, Path]="."):
    """Calls the function run_compute_electrostatics() to compute electrostatics potential of the PDB file.

    Args:
        pdb_filename (Union[str, Path]): Path to a PDB filename
        csv_filename (Union[str, Path]): Path to a CSV filename (output of run_compute_shell function)
        pqr_filename (Union[str, Path]): Optional path to a PQR filename
        out_dir (Union[str, Path]): Path to the output directory

    Returns:
        str: Path to a PDB file with coordinates of the shell particles and atomic electrostatic potential values in bfactor
    """


    # define access to useful APBS executables
    PATH_APBS = Path(os.getenv("APBS"))
    apbs = f"{str(PATH_APBS / 'bin' / 'apbs')}"
    inputgen = f"{str(PATH_APBS)}/share/apbs/tools/manip/inputgen.py"
    multivalue = f"{str(PATH_APBS)}/share/apbs/tools/bin/multivalue"

    # Create directory containing electrostatics calculation related files
    outdir = Path(out_dir)
    outdir.mkdir(exist_ok=True, parents=True)
    outfile_basename = str(outdir / Path(pdb_filename).stem)


    # compute pqr file if not given as input
    outfile_pqr = f"{outfile_basename}.pqr"
    if not pqr_filename:
        cmd_pdb2pqr = ["pdb2pqr30", "--ff", "CHARMM", "--whitespace", pdb_filename, outfile_pqr]
        logger.debug("Convert PDB to PQR format")
        status = subprocess.run(cmd_pdb2pqr, capture_output=True)
        if status != 0:
            logger.error(f"Error occured during the conversion of PDB to PQR, the process will stop.")
            exit(1)        
    else:
        logger.debug(f"Copying user-given PQR file as {outfile_pqr}")
        copyfile(src=pqr_filename, dst=outfile_pqr)


    # generate APBS input file
    outfile_in = f"{outfile_basename}.in"  # automaticall generated in the same path as the pqr file
    cmd_inputgen = ["python3", inputgen, "--method=auto", outfile_pqr]
    logger.debug(f"Running inputgen command: {' '.join(cmd_inputgen)}")
    status = subprocess.run(args=cmd_inputgen, capture_output=True)
    if status != 0:
        logger.error(f"Error occured during the inputgen command, the process will stop.")
        exit(1)


    # edit inputgen outfile
    logger.debug(f"Editing generated inputgen file: {outfile_in}")
    edit_inputgen(inputgen_file=outfile_in, pqrfile=outfile_pqr)


    # compute APBS electrostatics calculation
    outfile_pot = f"{outfile_pqr}.dx"  # default output name of APBS
    cmd_apbs = [apbs, outfile_in]
    logger.debug(f"Running APBS command: {cmd_apbs}")
    status = subprocess.run(args=cmd_apbs, capture_output=True)
    if status != 0:
        logger.error(f"Error occured during the APBS command, the process will stop.")
        exit(1)


    # cross the MSMS generated CSV file with the potential gridfile created by APBS
    outfile_multivalue = f"{outfile_basename}.mult"
    cmd_multivalue = [multivalue, csv_filename, outfile_pot, outfile_multivalue]
    logger.debug(f"Running multivalue command: {cmd_multivalue}")
    status = subprocess.run(args=cmd_multivalue, capture_output=True)
    if status != 0:
        logger.error(f"Error occured during the multivalue command, the process will stop.")
        exit(1)


    # convert .csv file into a PDB file
    outfile_pdb = f"{outfile_basename}_shell.pdb"
    logger.debug(f"Converting CSV to PDB file")
    run_multival_csv2pdb(ifile=outfile_multivalue, ofile=outfile_pdb)


    # remove io.mc generated by APBS
    io_apbs = Path.cwd() / "io.mc"
    logger.debug(f"Removing APBS derived log file: {io_apbs}")
    io_apbs.unlink(missing_ok=True)

    return outfile_pdb


def main():
    args = get_args()

    return run(
        pdb_filename=args.pdb,
        csv_filename=args.csv,
        pqr_filename=args.pqr,
        out_dir=args.out,
        out_subdir=args.subdir
    )

    
