#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Tuple, Union

from surfmap import PATH_MSMS, __COPYRIGHT_FULL__
from surfmap.lib.logs import get_logger
from surfmap.lib.parameters import Parameters
from surfmap.lib.utils import JunkFilePath
from surfmap.tools.SurfmapTools import run_particles_mapping
from surfmap.tools.compute_shell import run as run_compute_shell
from surfmap.tools.compute_electrostatics import run as run_compute_electrostatics


logger = get_logger(name=__name__)


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

    cmd = ["Rscript", params.coords_script, "-f", str(coords_file), "-s", str(params.cellsize), "-P", params.proj, "-o", params.outdir]
    
    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)
    
    if proc_status != 0:
        logger.error(f"Error occured during the computation of spherical coordinates of particles, the process will stop.")
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

    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)

    if proc_status != 0:
        logger.error(f"Error occured during smoothing of the raw matrix, the process will stop.")
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

    cmd = ["Rscript", params.map_script, "-i", matrix_file, scale_opt,  "-s", str(params.cellsize), "-p", params.pdb_id, "-P", params.proj, "-o", str(params.outdir), "--suffix", suffix]
    
    if params.coordstomap:
        cmd += ["-c", params.coordstomap]
    if reslist:
        cmd += ["-l", reslist]
    if params.png:
        cmd.append("--png")
        out_png = out_pdf.replace(".pdf", ".png")

    logger.debug(f"Running the command: {' '.join(cmd)}")
    proc_status = subprocess.call(cmd)

    if proc_status != 0:
        logger.error(f"Error occured during computing the map, the process will stop.")
        exit(1)

    return proc_status, out_png, out_pdf


def surfmap_from_pdb(params: Parameters, with_copyright: bool=True):
    """SURFMAP pipeline function to generate a 2D map from a PDB file. 

    Args:
        params (Parameters): an instance of Parameters object (-> surfmap.lib.parameters: class Parameters)
    """
    if with_copyright:
        print(__COPYRIGHT_FULL__)
    
    # create a trash for intermediary files to remove
    junk_optional = JunkFilePath(elements=[Path(params.outdir) / "shells", Path(params.outdir) / "tmp-elec"])
    shell = None

    listtomap = ["kyte_doolittle", "stickiness", "wimley_white", "circular_variance"] if params.ppttomap == "all" else [params.ppttomap]        
    for tomap in listtomap:

        logger.info(msg=f"Surface mapping of the {tomap} property".upper())

        step_index = 1
        if shell is None:
            outdir_shell = Path(params.outdir) / "shells"
            extra_radius = params.rad
            logger.info(msg=f"Step {step_index}: computing a shell around the protein surface")
            csv_coords, shell = run_compute_shell(pdb_filename=params.pdbarg, out_dir=outdir_shell, extra_radius=extra_radius)
        else:
            logger.info(msg=f"Step {step_index}: shell already exists. Skipping this step")


        if params.ppttomap == "electrostatics":
            step_index += 1
            outdir_elec = Path(params.outdir) / "tmp-elec"
            logger.info(msg=f"Step {step_index}: computing electrostatics potential")
            shell = run_compute_electrostatics(pdb_filename=params.pdbarg, csv_filename=csv_coords, pqr_filename=params.pqr, out_dir=outdir_elec)


        step_index += 1
        logger.info(msg=f"Step {step_index}: computing the property values and/or assign it to the shell particles")
        property = "bfactor" if tomap == "binding_sites" else tomap
        reslist, partlist_outfile = run_particles_mapping(shell=shell, pdb=params.pdbarg, tomap=property, outdir=params.outdir, res=params.resfile)
        junk_optional.add(element=[reslist, partlist_outfile])


        step_index += 1
        logger.info(msg=f"Step {step_index}: computing spherical coordinates of each shell particle and generate a raw matrix file")
        _, coordfile = compute_coords_list(params=params, coords_file=partlist_outfile, property=property)
        junk_optional.add(element=[coordfile, Path(coordfile).parent])

        
        step_index += 1
        logger.info(msg=f"Step {step_index}: smoothing the raw matrix file values")
        _, matrix_file, smoothed_matrix_file = compute_matrix(params=params, coords_file=coordfile, property=tomap)
        junk_optional.add(element=[matrix_file, Path(matrix_file).parent])


        step_index += 1
        logger.info(msg=f"Step {step_index}: computing the 2D map")
        _, png_filename, pdf_filename = compute_map(params=params, matrix_file=smoothed_matrix_file, property=tomap, reslist=reslist)


    # Deleting all intermediate files and directories
    if not params.keep:
        junk_optional.empty()
        
    # Creating log in output directory. Contains the parameters used to compute the maps.
    params.write_parameters(filename="parameters.log")
    print()


def surfmap_from_matrix(params: Parameters, with_copyright: bool=True):
    """SURFMAP function to generate a 2D map from a matrix file. 

    Args:
        params (Parameters): an instance of Parameters object (-> surfmap.lib.parameters: class Parameters)
    """
    if with_copyright:
        print(__COPYRIGHT_FULL__)

    if params.ppttomap == 'all':
        logger.error("Error: the property to map cannot be set to 'all' when computing a map from a matrix file.\n")
        exit()

    # create output dir and copy input matrix file in matrices/
    matrices_outdir = Path(params.outdir) / "matrices"
    matrices_outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(params.mat, matrices_outdir / Path(params.mat).name)

    matf = Path(params.outdir) / "matrices" / Path(params.mat).name
    
    logger.info(msg=f"Computing the 2D map")
    _, png_filename, pdf_filename = compute_map(params=params, matrix_file=matf, property=params.ppttomap)
    
    # Creating log in output directory. Contains the parameters used to compute the maps.
    params.write_parameters(filename="parameters.log")
 

