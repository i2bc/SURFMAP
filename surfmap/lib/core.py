#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
import shutil
import subprocess
from typing import Tuple, Union

from surfmap.lib.parameters import Parameters
from surfmap.lib.utils import JunkFilePath


def compute_shell(params: Parameters) -> Tuple[int, str]:
    """Compute the shell of a given PDB structure.

    Two possible options: simple shell or shell with electrostatic potential computed by APBS.

    Args:
        params (Parameters): Object with a set of useful parameters:
            - params.ppttomap  # property to map
            - params.shell_script  # path to the script compute_shell.sh
            - params.pdbarg  # pdb filename
            - rad: float  # radius used for shell computation
            - params.outdir  # path to the output directory

    Returns:
        tuple[int, str]: return the status value of the process and the output filename
    """
    elecval = "1" if params.ppttomap == "electrostatics" else "0"

    # Create shell. (two options: simple shell or shell with elec potential computed by APBS.)
    proc_status = subprocess.call(["bash", params.shell_script, "-p", params.pdbarg, "-e", elecval, "-r", str(params.rad), "-o", params.outdir])
    outfile = str(Path(params.outdir) / "shells" / (Path(params.pdbarg).stem + "_shell.pdb"))

    if proc_status != 0:
        print("The script compute_shell.sh failed. It can be because MSMS could not compute the Connolly surface of the protein.\nYou need to provide a pdb file that MSMS can handle.\nPlease check that there is no problem with your input pdb file.\nExiting now.")
        exit()
    if not Path(outfile).exists():
        print("Shell not found. This is probably because MSMS did not manage to compute the surface of the protein.\nExiting now.")
        exit()

    return proc_status, outfile


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
    out_coords_reslist = None  # output file with spherical coordinates of specific residues
    out_coords_all = Path(params.outdir) / f"{Path(params.pdbarg).stem}_{property}_partlist.out"  # output file with spherical coordinates of each residue
    out_pdb_cv = Path(params.curdir) / f"{Path(params.pdbarg).stem}_CV.pdb" 

    if params.resfile:
        cmd_surfmap = [params.surftool_script, "-pdb", params.pdbarg, "-shell", shell_filename, "-tomap", property, "-res", params.resfile, "-d", params.outdir]
        out_coords_reslist = str(Path(params.outdir) / Path(params.resfile).stem + "_sph_coords.out")
    else:
        cmd_surfmap = [params.surftool_script, "-pdb", params.pdbarg, "-shell", shell_filename, "-tomap", property, "-d", params.outdir]

    proc_status = subprocess.call(cmd_surfmap)
    if out_pdb_cv.exists():
        # out_pdb_cv.unlink()
        try:
            shutil.move(str(out_pdb_cv), str(params.outdir))
        except OSError:
            pass
    
    if proc_status != 0:
        print("\nMapping of the residues given in input failed. This can be due to a malformed residue file or a mistake in the reisdue numbering, type or chain.\n\nFormat should be the following:\nchain residuenumber residuetype\nexample: A\t5 LEU\n\nExiting now.")
        exit()

    return proc_status, out_coords_reslist, out_coords_all


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
    out_file = str(Path(params.outdir) / "coord_lists" / f"{Path(params.pdbarg).stem}_{property}_coord_list.txt")

    cmdlist = ["Rscript", params.coords_script, "-f", coords_file, "-s", str(params.cellsize), "-P", params.proj]
    proc_status = subprocess.call(cmdlist)

    return proc_status, out_file


def compute_matrix(params: Parameters, coords_file: str, property: str) -> Tuple[int, str, str]:
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
    out_matrix_smoothed = str(Path(params.outdir) / "smoothed_matrices" / f"{Path(params.pdbarg).stem}_{named_property}_smoothed_matrix.txt")
    out_matrix = str(Path(params.outdir) / "matrices" / f"{Path(params.pdbarg).stem}_{named_property}_matrix.txt")

    cmd = ["Rscript", params.matrix_script, "-i", coords_file, "-s", str(params.cellsize), "--discrete", "-P", params.proj]
    if params.nosmooth and property != "binding_sites":
        cmd[-3] = "--nosmooth"
    elif not params.nosmooth and property != "binding_sites":
        del cmd[-3]

    proc_status = subprocess.call(cmd)

    return proc_status, out_matrix, out_matrix_smoothed        


def compute_map(params: Parameters, matrix_file: str, property: str, reslist: str=None) -> int:
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
    out_pdf = str(Path(params.outdir) / "maps" / f"{params.pdb_id}_{property}_map.pdf")
    out_png = None

    if property == "binding_sites":
        scale_opt = "--discrete"
    elif property == "circular_variance_atom":
        scale_opt = "--circular_variance"
    else:
        scale_opt = "--" + property

    cmdmap = ["Rscript", params.map_script, "-i", matrix_file, scale_opt,  "-s", str(params.cellsize), "-p", params.pdb_id, "-P", params.proj]
    
    if params.coordstomap:
        cmdmap += ["-c", params.coordstomap]
    if reslist:
        cmdmap += ["-l", reslist]
    if params.png:
        cmdmap.append("--png")
        out_png = out_pdf.replace(".pdf", ".png")

    proc_status = subprocess.call(cmdmap)

    return proc_status, out_png, out_pdf



def surfmap_from_pdb(params: Parameters):
    # create a junk for intermediary files to remove
    junk_optional = JunkFilePath(elements=[Path(params.outdir) / "shells", Path(params.outdir) / "tmp-elec"])


    # Step 1: Generation of shell around protein surface.
    _, shell = compute_shell(params=params)

    listtomap = ["kyte_doolittle", "stickiness", "wimley_white", "circular_variance"] if params.ppttomap == "all" else [params.ppttomap]        
    for tomap in listtomap:        
        print("Surface property mapping:", tomap)

        # Step 2: compute and/or associate the property value of interest (electrostatics, hydrophobicity, stickiness...) to atoms/residues
        property = "bfactor" if tomap == "binding_sites" else tomap
        _, reslist, mapfile = compute_spherical_coords(params=params, shell_filename=shell, property=property)
        junk_optional.add(element=[reslist, mapfile])

        # Part 3: compute phi theta list and generete a raw matrix file
        _, coordfile = compute_coords_list(params=params, coords_file=mapfile, property=property)
        junk_optional.add(element=[coordfile, Path(coordfile).parent])

        # Step 4: smooth the matrix
        _, matrix_file, smoothed_matrix_file = compute_matrix(params=params, coords_file=coordfile, property=tomap)
        junk_optional.add(element=[matrix_file, Path(matrix_file).parent])

        # Step 5: computing map
        _, png_filename, pdf_filename = compute_map(params=params, matrix_file=smoothed_matrix_file, property=tomap, reslist=reslist)


    # Deleting all intermediate files and directories
    if not params.keep:
        junk_optional.empty()
        
    # Creating log in output directory. Contains the parameters used to compute the maps.
    params.write_parameters(filename="parameters.log")


def surfmap_from_matrix(params: Parameters):
    """Compute SURFMAP from a matrix file and a property to map

    Args:
        params (Parameters): Instance with a set of useful parameters. (-> surfmap.lib.parameters: class Parameters)
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
 

