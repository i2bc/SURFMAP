
import importlib.resources as pkg_resources
from pathlib import Path
import os

with pkg_resources.path("surfmap.scripts", 'compute_shell.sh') as resource:
    ROOT = Path(resource).parent.parent.parent.resolve()

PATH_R_SCRIPTS = ROOT / "surfmap" / "scripts"
PATH_MSMS = ROOT / "surfmap" / "utils" / "MSMS"
PATH_TO_EXAMPLES = ROOT / "surfmap"/ "examples"


__VERSION__ = '2.0.0'
SURFMAP_DOCKER_VERSION = os.getenv('SURFMAP_DOCKER_VERSION', __VERSION__)

SEP = f"{'-' * 80}"

__COPYRIGHT_NOTICE__ = """
SURFMAP: Projection of protein surface features on 2D map
Copyright (C) 2021  H. Schweke
Version: {}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
{}""".format(__VERSION__, SEP)


__COPYRIGHT_MSMS__ = """
SURFMAP relies on the use of MSMS (version 2.6.1) - Copyright (c) M. F. Sanner, 1994.

As a corollary of MSMS's license terms and conditions, by using SURFMAP,
you agree to acknowledge the use of the MSMS software that results in any 
published work, including scientific papers, films and videotapes by 
citing the following reference:

    Sanner, M.F., Spehner, J.-C., and Olson, A.J. (1996) Reduced surface:
    an efficient way to compute molecular surfaces. Biopolymers, Vol. 38.,
    (3), 305-320.

MSMS is free for academic use. For commercial use please contact
M. F. Sanner at sanner@scripps.edu.
{}
If you use SURFMAP for your research, please cite the following
papers:
    - Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
      Simon Gosset, Anne Lopes (2021) SURFMAP: a software for mapping
      in two dimensions protein surface features. bioRxiv
      2021.10.15.464543; doi: https://doi.org/10.1101/2021.10.15.464543

    - Sanner, M.F., Spehner, J.-C., and Olson, A.J. (1996) Reduced 
      surface: an efficient way to compute molecular surfaces. 
      Biopolymers, Vol. 38., (3), 305-320.
{}
""".format(SEP, SEP)

__COPYRIGHT_FULL__ = f"{__COPYRIGHT_NOTICE__}{__COPYRIGHT_MSMS__}"