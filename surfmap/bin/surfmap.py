#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import surfmap
from surfmap import __COPYRIGHT_FULL__
from surfmap.lib.core import surfmap_from_pdb, surfmap_from_matrix
from surfmap.lib.logs import set_root_logger
from surfmap.lib.docker import DockerCLI
from surfmap.lib.parameters import Parameters, get_args


def surfmap_local(params: Parameters):
    """Execute SURFMAP from a local install

    Args:
        params (Parameters): Set of useful parameters.
    """
    if params.pdbarg:
        surfmap_from_pdb(params=params)
    elif params.mat:
        surfmap_from_matrix(params=params)


def surfmap_container(params: Parameters):
    """Execute SURFMAP from a docker container.

    This requires to convert the command arguments in order to:
        1. map required input/output mounting points between local and container file systems
        2. run the container as a subprocess

    Args:
        params (Parameters): Set of useful parameters.
    """
    if params.pdbarg:
        input_args = ['-pdb']
        if params.ppttomap == "electrostatics" and params.pqr:
            input_args.append('-pqr')
        if params.resfile:
            input_args.append('-res')
    elif params.mat:
        input_args = ['-mat']


    cli = DockerCLI(
        input_args=input_args,
        output_args="-d",
        out_dirname=params.outdir,
        input_dir="/home/surfmap/input",
        output_dir="/home/surfmap/output",
        args=params.args
    )

    cli.show()
    cli.run()


def main():
    params = Parameters(args=get_args())
    set_root_logger(outpath=params.outdir, level=params.verbose)

    if params.docker:
        surfmap_container(params=params)
    else:
        surfmap_local(params=params)


if __name__ == "__main__":
    main()