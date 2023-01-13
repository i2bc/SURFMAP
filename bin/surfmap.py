#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from surfmap.lib.parameters import Parameters
from surfmap.lib.computing import surfmap_from_pdb, surfmap_from_matrix


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
       
    
def main():
    params = Parameters()
    show_copyrights()

    if params.pdbarg:
        surfmap_from_pdb(params=params)
    elif params.mat:
        surfmap_from_matrix(params=params)


if __name__ == "__main__":
    main()

