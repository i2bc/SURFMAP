#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess

from surfmap.lib.parameters import Parameters
from surfmap.lib.core import surfmap_from_pdb, surfmap_from_matrix
from surfmap.lib.docker import DockerCLI


def show_copyrights():
    msg = """
SURFMAP:    Projection of protein surface features on 2D map
Authors:    Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
            Simon Gosset, Anne Lopes
Version:    1.5
Copyright (c) 2021, H. Schweke


SURFMAP relies on the use of MSMS (Copyright (c) M. F. Sanner, 1994).
MSMS is free for academic use. For commercial use please contact
M. F. Sanner at sanner@scripps.edu.

If you use SURFMAP for your research, please cite the following
papers:
    - Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
      Simon Gosset, Anne Lopes (2021) SURFMAP: a software for mapping
      in two dimensions protein surface features. bioRxiv
      2021.10.15.464543; doi: https://doi.org/10.1101/2021.10.15.464543

    - Sanner, M.F., Spehner, J.-C., and Olson, A.J. (1996) Reduced 
      surface: an efficient way to compute molecular surfaces. 
      Biopolymers, Vol. 38., (3), 305-320.

SURFMAP can also optionnaly map electrostatic potential through the use
of the APBS software. So if you use SURFMAP with the option 
'-tomap electrostatics' and find it useful to your research please cite
one or more of the papers listed in the following webpage:
https://apbs.readthedocs.io/en/latest/supporting.html
---------------------------------------------------------------------------

    """
    print(msg)



def surfmap_local(params: Parameters):
    """Execute SURFMAP from a local install

    Args:
        params (Parameters): Set of useful parameters.
    """
    show_copyrights()

    if params.pdbarg:
        surfmap_from_pdb(params=params)
    elif params.mat:
        surfmap_from_matrix(params=params)


def surfmap_container(params: Parameters):
    """Execute SURFMAP from a docker container.

    This requires to convert the command arguments in order to:
        1. map required input/output mount point between local and container file systems
        2. run the container as sa subprocess

    Args:
        params (Parameters): Set of useful parameters.
    """
    cli = DockerCLI(
        input_args=["-pdb"],
        output_args="-d",
        out_dirname=params.outdir,
        input_dir="/home/surfmap/input",
        output_dir="/home/surfmap/output",
        args=params.args
    )

    cli.show()
    cli.run()


def main():
    params = Parameters()

    if params.docker:
        surfmap_container(params=params)
    else:
        surfmap_local(params=params)


if __name__ == "__main__":
    main()