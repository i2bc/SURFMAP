#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import shutil
import subprocess
from typing import Tuple, Union

from surfmap import PATH_MSMS
from surfmap.lib.custom_logging import Log
from surfmap.lib.parameters import Parameters
from surfmap.lib.utils import JunkFilePath
from surfmap.tools.SurfmapTools import run_spherical_coords
from surfmap.tools.compute_shell import run as run_compute_shell
from surfmap.tools.compute_electrostatics import run as run_compute_electrostatics


logger = Log()


def compute_shell(params: Parameters) -> Tuple[int, str]:
    """Compute the shell of a given PDB structure. If '-tomap electrostatics' property is requested, 
    then electrostatics calculation is done during this step.

    Two possible options: simple shell or shell with electrostatic potential computed by APBS.

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.ppttomap  # property to map
            - params.shell_script  # path to the script compute_shell.sh
            - params.pdbarg  # pdb filename
            - params.rad: float  # radius used for shell computation
            - params.pqr: str  # path to a pqr file (optionally used with electrostatics only)
            - params.outdir  # path to the output directory

    Returns:
        tuple[int, str]: return the status value of the process and the output filename
    """

    # compute shell surface - returns a PDB file with values set to 0.000 in bfactor 
    outfile_csv, outfile_pdb = run_compute_shell(
            pdb_filename=params.pdbarg,
            out_dir=params.outdir,
            out_subdir="shells",
            extra_radius=params.rad
    )

    # compute electrostatics of a protein surface - returns a PDB file with electrostatics values in bfactor 
    if params.ppttomap == "electrostatics":
        print(" - computing electrostics")
        outfile_pdb = run_compute_electrostatics(
            pdb_filename=params.pdbarg,
            csv_filename=outfile_csv,
            out_dir=params.outdir,
            out_subdir="tmp-elec"
        )

    return outfile_pdb


def compute_shell_old(params: Parameters) -> Tuple[int, str]:
    """Compute the shell of a given PDB structure.

    Two possible options: simple shell or shell with electrostatic potential computed by APBS.

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.ppttomap  # property to map
            - params.shell_script  # path to the script compute_shell.sh
            - params.pdbarg  # pdb filename
            - rad: float  # radius used for shell computation
            - pqr: str  # path to a pqr file (optionally used with electrostatics only)
            - params.outdir  # path to the output directory

    Returns:
        tuple[int, str]: return the status value of the process and the output filename
    """
    outfile_shell = str(Path(params.outdir) / "shells" / (Path(params.pdbarg).stem + "_shell.pdb"))  # an output file generated by params.shell_script
    elecval = "1" if params.ppttomap == "electrostatics" else "0"

    # Create shell. (two options: simple shell or shell with elec potential computed by APBS.)
    cmd_shell = ["bash", params.shell_script, "-p", params.pdbarg, "-e", elecval, "-r", str(params.rad), "-o", params.outdir]

    if params.pqr:
        cmd_shell += ["-q", params.pqr]

    proc_status = subprocess.call(cmd_shell)

    if proc_status != 0:
        print("The script compute_shell.sh failed. It can be because MSMS could not compute the Connolly surface of the protein.\nYou need to provide a pdb file that MSMS can handle.\nPlease check that there is no problem with your input pdb file.\nExiting now.")
        exit()
    if not Path(outfile_shell).exists():
        print("Shell not found. This is probably because MSMS did not manage to compute the surface of the protein. Exiting now.\n")
        exit()

    return proc_status, outfile_shell


def compute_spherical_coords(params: Parameters, shell_filename: str, property: str) -> Tuple[int, str, str]:
    """Convert cartesian coordinates of each particule into spherical coordinates and 
    associates the value of interest (electrostatics, hydrophobicity, stickiness...)

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.resfile  # residue filename to map, if given
            - params.ppttomap  # property to map
            - params.surftool_script  # path to the executable _surfmap_tool (-> surfmap.tools.SurfmapTools)
            - params.pdbarg  # pdb filename
            - params.outdir  # path to the output directory
        
        shell_filename (str): Path to the filename of a computed shell
        property (str): property to map

    Returns:
        tuple[int, str, str]: return the status value of the process, the output filename of reslist, if any, and the output filename of partlist
    """
    outfile_res_to_map, partlist_outfile = run_spherical_coords(shell=shell_filename, pdb=params.pdbarg, tomap=property, outdir=params.outdir, res=params.resfile)

    return outfile_res_to_map, partlist_outfile


def compute_coords_list(params: Parameters, coords_file: str, property: str) -> Tuple[int, str]:
    """Compute phi theta list

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.coords_script  # path to the script computeCoordsList.r
            - params.pdbarg  # pdb filename
            - params.outdir  # path to the output directory
            - params.cellsize  # unit size of a grid cell
            - params.proj  # name of the projection type

        coords_file (str): # path to filename with spherical coordinates of each residue
        property (str): property to map

    Returns:
        tuple[int, str]: return the status value of the process and the output filename
    """
    out_file = str(Path(params.outdir) / "coord_lists" / f"{Path(params.pdbarg).stem}_{property}_coord_list.txt")  # an output file generated by params.coords_script

    cmdlist = ["Rscript", params.coords_script, "-f", coords_file, "-s", str(params.cellsize), "-P", params.proj, "-o", params.outdir]
    proc_status = subprocess.call(cmdlist)

    if proc_status != 0:
        print(f"Error occured in step 3, the process will stop.")
        exit(1)

    return proc_status, out_file


def compute_matrix(params: Parameters, coords_file: str, property: str, suffix="_coord_list.txt") -> Tuple[int, str, str]:
    """Compute matrix files

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.matrix_script  # path to the script computeMatrices.r
            - params.pdbarg  # pdb filename
            - params.outdir  # path to the output directory
            - params.cellsize  # unit size of a grid cell
            - params.proj  # name of the projection type

        coords_file (str): # path to filename with coords_list (output of compute_coords_list())
        property (str): property to map

    Returns:
        tuple[int, str, str]: return the status value of the process, the output filenames of the matrix file and the smoothed matrix file.
    """
    named_property = "bfactor" if property == "binding_sites" else property
    out_matrix_smoothed = str(Path(params.outdir) / "smoothed_matrices" / f"{Path(params.pdbarg).stem}_{named_property}_smoothed_matrix.txt")  # an output file generated by params.matrix_script
    out_matrix = str(Path(params.outdir) / "matrices" / f"{Path(params.pdbarg).stem}_{named_property}_matrix.txt")  # an output file generated by params.matrix_script

    cmd = ["Rscript", params.matrix_script, "-i", coords_file, "-s", str(params.cellsize), "-P", str(params.proj), "-o", str(params.outdir), "--suffix", suffix, "--discrete"]

    if params.nosmooth and property != "binding_sites":
        cmd[-1] = "--nosmooth"
    elif not params.nosmooth and property != "binding_sites":
        del cmd[-1]

    proc_status = subprocess.call(cmd)

    if proc_status != 0:
        print(f"Error occured in step 4, the process will stop.")
        exit(1)

    return proc_status, out_matrix, out_matrix_smoothed        


def compute_map(params: Parameters, matrix_file: str, property: str, reslist: str=None, suffix="_smoothed_matrix.txt") -> int:
    """Compute map from a matrix file

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.map_script  # path to the script computeMaps.r
            - params.pdb_id  # PDB stem name (e.g. '1g3n')
            - params.coordstomap  # None
            - params.outdir  # path to the output directory
            - params.cellsize  # unit size of a grid cell
            - params.proj  # name of the projection type
        matrix_file (str): Path to a smoothed matrix file.
        property (str): property to map
        reslist (str): Path to the file with spherical coordinates of specific residues

    Returns:
        int: status value of the process, png output file and pdf output file 
    """
    out_pdf = str(Path(params.outdir) / "maps" / f"{params.pdb_id}_{property}_map.pdf")  # an output file generated by params.map_script
    out_png = None

    if property == "binding_sites":
        scale_opt = "--discrete"
    elif property == "circular_variance_atom":
        scale_opt = "--circular_variance"
    else:
        scale_opt = "--" + property

    cmdmap = ["Rscript", params.map_script, "-i", matrix_file, scale_opt,  "-s", str(params.cellsize), "-p", params.pdb_id, "-P", params.proj, "-o", str(params.outdir), "--suffix", suffix]
    
    if params.coordstomap:
        cmdmap += ["-c", params.coordstomap]
    if reslist:
        cmdmap += ["-l", reslist]
    if params.png:
        cmdmap.append("--png")
        print(cmdmap)
        out_png = out_pdf.replace(".pdf", ".png")

    proc_status = subprocess.call(cmdmap)

    if proc_status != 0:
        print(f"Error occured in step 5, the process will stop.")
        exit(1)

    return proc_status, out_png, out_pdf



def surfmap_from_pdb(params: Parameters):
    """SURFMAP pipeline function to generate a 2D map from a PDB file. 

    Args:
        params (Parameters): an instance of Parameters object (-> surfmap.lib.parameters: class Parameters)
    """

    # create a trash for intermediary files to remove
    junk_optional = JunkFilePath(elements=[Path(params.outdir) / "shells", Path(params.outdir) / "tmp-elec"])
    shell = None

    listtomap = ["kyte_doolittle", "stickiness", "wimley_white", "circular_variance"] if params.ppttomap == "all" else [params.ppttomap]        
    for tomap in listtomap:

        logger.section(message=f"Surface mapping of the {tomap} property".upper())

        if shell is None:
            logger.info(message="- Step 1: computing a shell around the protein surface")
            shell = compute_shell(params=params)
        else:
            logger.info(message="- Step 1: shell already exists. Skipping this step")

        logger.info(message="- Step 2: computing the property values and/or assign it to the shell particles")
        property = "bfactor" if tomap == "binding_sites" else tomap
        reslist, partlist_outfile = compute_spherical_coords(params=params, shell_filename=shell, property=property)
        junk_optional.add(element=[reslist, partlist_outfile])


        logger.info(message="- Step 3: computing spherical coordinates of each shell particle and generate a raw matrix file")
        _, coordfile = compute_coords_list(params=params, coords_file=partlist_outfile, property=property)
        junk_optional.add(element=[coordfile, Path(coordfile).parent])


        logger.info(message="- Step 4: Smoothing the raw matrix file values")
        status, matrix_file, smoothed_matrix_file = compute_matrix(params=params, coords_file=coordfile, property=tomap)
        junk_optional.add(element=[matrix_file, Path(matrix_file).parent])


        logger.info(message="- Step 5: Computing the 2D map")
        _, png_filename, pdf_filename = compute_map(params=params, matrix_file=smoothed_matrix_file, property=tomap, reslist=reslist)


    # Deleting all intermediate files and directories
    if not params.keep:
        junk_optional.empty()
        
    # Creating log in output directory. Contains the parameters used to compute the maps.
    params.write_parameters(filename="parameters.log")
    print()


def surfmap_from_matrix(params: Parameters):
    """SURFMAP pipeline function to generate a 2D map from a matrix file. 

    Args:
        params (Parameters): an instance of Parameters object (-> surfmap.lib.parameters: class Parameters)
    """

    if params.ppttomap == 'all':
        print("\nError: the property to map cannot be set to 'all' when computing a map from a matrix file.\n".format(params.mat))
        exit()

    # create output dir and copy input matrix file in matrices/
    matrices_outdir = Path(params.outdir) / "matrices"
    matrices_outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(params.mat, matrices_outdir / Path(params.mat).name)

    matname =  Path(params.mat).stem
    matf = Path(params.outdir) / "matrices" / Path(params.mat).name
    
    # Step 5: computing map        
    _, png_filename, pdf_filename = compute_map(params=params, matrix_file=matf, property=params.ppttomap)
    
    # Creating log in output directory. Contains the parameters used to compute the maps.
    params.write_parameters(filename="parameters.log")
 

