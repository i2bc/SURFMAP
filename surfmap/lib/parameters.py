import argparse
from pathlib import Path
from shutil import which
import sys
from typing import Any, Union

from surfmap import PATH_R_SCRIPTS, __COPYRIGHT_NOTICE__


def get_args():
    parser = argparse.ArgumentParser()
    group_mutually_exclusive = parser.add_mutually_exclusive_group(required=True)
    
    group_mutually_exclusive.add_argument(
        "-pdb",
        help="Input pdb file (path + file name)"
    )

    group_mutually_exclusive.add_argument(
        "-mat",
        type=str,
        help="Input matrix. If the user gives an imput matrix, SURFMAP will directly compute a map from it."
    )

    group_mutually_exclusive.add_argument(
        "-v", "--version",
        action="store_true",
        help="Print the current used version of SURFMAP."
    )

    parser.add_argument(
        "-tomap",
        type=str,
        required=True,
        choices=set(("all", "stickiness", "kyte_doolittle", "wimley_white", "electrostatics", "circular_variance", "circular_variance_atom", "bfactor", "binding_sites")),
        help="Choice of the scale. Argument must be one of the following: stickiness; kyte_doolittle; wimley_white; electrostatics; circular_variance; bfactor; binding_sites; all"
    )

    parser.add_argument(
        "-proj",
        type=str,
        required=False,
        default="flamsteed",
        choices=set(("flamsteed", "mollweide", "lambert")),
        help="Choice of the projection. Argument must be one of the following: flamsteed; mollweide; lambert (default: flamsteed)"
    )


    parser.add_argument(
        "-coords",
        default=None,
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        "-res",
        type=str,
        default='',
        help="File containing a list of residues to map on the projection. Format must be the following: col 1 = chain id; col 2 = res number; col 3 = res type"
    )

    parser.add_argument(
        "-rad",
        required=False,
        type=float,
        default=3.0,
        help="Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the smoother the surface (default: 3.0)"
    )

    parser.add_argument(
        "-d",
        required=False,
        type=str,
        default="output_SURFMAP_{}_{}",
        help = "Output directory where all files will be written (default: './output_SURFMAP_$pdb_$tomap' where $pdb and $tomap are the inputs given to -pdb and -tomap arguments, respectiveley)"
    )

    parser.add_argument(
        "-s",
        required=False,
        type=int,
        default=5,
        help="Size of a grid cell. Necessary that 180%%cellsize == 0 (default: 5.0)"
    )

    parser.add_argument(
        "--nosmooth",
        action="store_true",
        help="If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)"
    )

    parser.add_argument(
        "--png",
        action="store_true", 
        help="If chosen, a map in png format is computed (default: only pdf format is generated)"
    )

    parser.add_argument(
        "--keep",
        action="store_true",
        help="If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)"
    )

    parser.add_argument(
        "--docker",
        action="store_true",
        help="If chosen, SURFMAP will be run on a docker container (requires docker installed)."
    )

    parser.add_argument(
        "--pqr",
        required=False,
        type=str,
        default=None,
        help="Path to a PQR file used for electrostatics calculation. Option only available if '-tomap electrosatics' is requested. Defaults to None."
    )

    parser.add_argument(
        "-verbose",
        required=False,
        type=int,
        default=1,
        help="Verbose level of the console log. 0 for silence, 1 for info level, 2 for debug level. Defaults to 1."
    )
    
    if True in [ x in ["-v", "--version"] for x in sys.argv[1:] ]:
        print(f"{__COPYRIGHT_NOTICE__}")
        exit()

    return parser.parse_args()


class Parameters:
    """Class handler for parameters of surfmap 
    
    Attributes are:
    - curdir: Union[str, Path]  # current working directory
    - surftool_script: str  # path to binary _surfmap_tool 
    - shell_script: str  # path to compute_shell.sh
    - coords_script: str  # path to computeCoordList.sh
    - matrix_script: str  # path to computeMatrices.sh
    - map_script: str  # path to computeMaps.sh
    - mat: str=None  # path to a given matrix file
    - pdbarg: str = None  # path to a given PDB file
    - pdb_id: str  # PDB stem name (e.g. '1g3n')
    - pdbname: str  # name of PDB file with no path (e.g. '1g3n.pdb')
    - proj: str  # name of the projection type
    - ppttomap: # name of the property
    - resfile: str  # residue filename to map, if given
    - rad: float  # radius used for shell computation
    - cellsize: str  # unit size of a grid cell
    - coordstomap: Any = None  #
    - nosmooth: bool  # True to have map not smoothed (for discrete values only)
    - png: bool  # True to generate a PNG file of the map in addition to the usual PDF
    - keep: bool  # True to keep intermediary files that are usually removed
    - docker: bool  # True to run SURFMAP on a docker container
    - outdir: Union[str, Path]  # path to the output directory
    - pqr: str=None  # path to a PQR file used for electrostatics calculation. Defaults to None
    - verbose: int=2  # Verbose level of the console log. 0 for silence, 1 for debug level, 2 for info level. Defaults to 2.

    """
    REQUIREMENTS = ["R", "awk", "apbs"]

    PROJECTION_MAP = {'flamsteed': 'sinusoidal', 'mollweide': 'mollweide', 'lambert': 'lambert'}
    PROPERTIES = {"all": ["kyte_doolittle", "stickiness", "wimley_white", "circular_variance"]}

    VERBOSE_MAP = {
        0: 100,
        1: 20,
        2: 10,
        3: 30,
        4: 40,
        5: 50
    }
    
    def __init__(self, args: Union[argparse.Namespace, object], path_to_scripts: Union[str, Path]=PATH_R_SCRIPTS) -> None:
        self.args = args

        if not args.docker:
            self._check_surfmap_requirements()
        else:
            self._check_docker_install()

        # define pdb relative variables
        if args.pdb:
            self.mat: str = None
            self.pdbarg: str = args.pdb
            self._check_pdbarg()
            self.pdb_id: str = Path(self.pdbarg).stem if not args.mat else Path(args.mat).stem
            self.pdbname: str = Path(self.pdbarg).name

        # define matrice if given in input
        elif args.mat:
            self.pdbarg: str = None
            self.pdb_id: str = Path(args.mat).stem
            self.mat: str = args.mat
            self.pdbname: str = Path(self.mat).name

        # set useful paths
        self.curdir: Union[str, Path] = Path.cwd()
        self.surftool_script: str = "_surfmap_tool"
        self.shell_script: str = str(Path(path_to_scripts) / "compute_shell.sh")
        self.coords_script: str  = str(Path(path_to_scripts) / "computeCoordList.R")
        self.matrix_script: str = str(Path(path_to_scripts) / "computeMatrices.R")
        self.map_script: str = str(Path(path_to_scripts) / "computeMaps.R")

        # define projection type
        self.proj: str = self.PROJECTION_MAP[args.proj]

        # define property to map
        self.ppttomap: str = args.tomap

        # define optional pqr file (optionally used for electrostatics)
        self.pqr: str = args.pqr

        # define residues to map, if any
        self.resfile: str = args.res
        if self.resfile and not Path(self.resfile).exists():
            print("The residue file could not be found (arg -res). It seems that this file does not exist.\nThis could be due to a mistake in the path to the file, for example.\nExiting now.")
            exit()

        # define radius in angström added to usual atomic radius (used for calculation solvent excluded surface)
        self.rad: float = args.rad
        self._check_rad()
        
        # define cellsize of the map
        self.cellsize: int = args.s
        self._check_cellsize()

        # define coords to map; not used
        self.coordstomap: Any = args.coords

        self.nosmooth: bool = args.nosmooth
        self.png: bool = args.png
        self.keep: bool = args.keep
        self.docker: bool = args.docker
        
        # define and create output directory if not exists
        self.DEFAULT_OUTDIR_BASENAME = "output_SURFMAP_{}_{}"
        self._set_outdir(args=args)

        # set verbose
        if args.verbose in self.VERBOSE_MAP:
            self.verbose = self.VERBOSE_MAP[args.verbose]
        else:
            print("Warning: verbose level must be either 0, 1, or 2. Use of default verbose level (2)")
            self.verbose = self.VERBOSE_MAP[2]

    def _check_surfmap_requirements(self):
        """Check if requirements are satisfied (will exit if not).

        Only consider APBS if 'electrostatics' is asked as a -tomap option argument.
        """
        exe_not_found = []
        for requirement in self.REQUIREMENTS:
            if requirement == "apbs" and self.args.tomap != "electrostatics":
                continue
            if not which(requirement):
                exe_not_found.append(requirement)
            
        if exe_not_found:
            print(f"\nError: the following executable(s) required to run SURFMAP has(ve) not been detected on your system: {', '.join(exe_not_found)}.\nPlease install the missing executable(s).\n")
            print(f"If you think that an executable is present on your system but is detected as missing, please make sure to make it accessible no matter the current directory. It can be done by setting your search path (export PATH=$PATH:/...)\n")
            exit()

    
    def _check_docker_install(self):
        if not which("docker"):
            print(f"Error: docker has not been detected on your system. Please install it to use a pre-built image of SURFMAP.\n")


    def _check_mutually_exclusive_args(self, args):
        if not args.pdb and not args.mat:
            print("You need to provide either a pdb file or a matrix file in input.\nexiting now.")
            exit()
        elif args.pdb and args.mat:
            print("You need to provide either a pdb file or a matrix file in input. you cannot provide both.\nexiting now.")
            exit()

    def _check_pdbarg(self):
        if not Path(self.pdbarg).exists():
            print("Error: the pdb file {} does not exist".format(self.pdbarg))
            exit()
        if Path(self.pdbarg).suffix != ".pdb":
            print("The file provided in input seems to not be a pdb file. Please provide a file with a '.pdb' extension.")
            exit()

    def _check_rad(self):
        if not 10.0 >= self.rad >= 0.5:
            print("The radius given in input (arg -rad) is either too large or too small. Please provide a number between 0.5 and 10.\nExiting now.")
            exit()

    def _check_cellsize(self):
        if 180 % self.cellsize != 0:
            print("cellsize (arg -s) is not a multiplier of 180.\nThe cell size needs to be a multiplier of 180 in order to compute the grid. \nExiting now.")
            exit()
        elif self.cellsize > 20:
            print("the cellsize is too large (arg -s).\nThe maximum cell size authorized is 20, which represents a 9*18 squares grid. \nExiting now.")
            exit()

    def _set_outdir(self, args):
        self._is_outdir_default = True

        if args.d != self.DEFAULT_OUTDIR_BASENAME:
            self.outdir: Union[str, Path] = args.d
            self._is_outdir_default = False
        else:
            names = [self.pdb_id, 'all_properties'] if self.ppttomap == "all" else [self.pdb_id, self.ppttomap]
            self.outdir: Union[str, Path] = args.d.format(*names)

        Path(self.outdir).mkdir(parents=True, exist_ok=True)
                
    def get_log_parameters(self):
        return """
Parameters used to compute the maps:
- PDB file (-pdb): {}
- Matrix file (-mat): {}
- Projection type (-proj): {}
- PQR file: {}
- Name of the property mapped (-ppttomap): {} 
- Filename of residues to map (-resfile): {}
- MSMS radius used for shell computation (-rad): {}
- Unit size of the grid cell (-s): {}
- Grid resolution: {}
- Map not smoothed (--nosmooth): {}
- Generate png (--png): {}
- Keep intermediary files (--keep): {}
- Output directory (-outdir): {}
    """.format(
        self.pdbarg if self.pdbarg else "None",
        self.mat if self.mat else "None",
        self.proj,
        self.pqr if self.pqr else "None",
        self.ppttomap if self.ppttomap != "all" else ", ".join(self.PROPERTIES["all"]),
        self.resfile if self.resfile else "None",
        self.rad,
        self.cellsize,
        f"{int(360 / self.cellsize)}*{int(180 / self.cellsize)}",
        self.nosmooth,
        self.png,
        self.keep,
        self.outdir
    )

    def write_parameters(self, filename: str, outdir: Union[str, Path]=None):
        outdir = Path(self.outdir) if not outdir else Path(outdir)
        out_filename = outdir / filename

        with open(out_filename, "w+") as out_file:
            out_file.write(self.get_log_parameters())
