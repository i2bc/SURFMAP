<div align="center">

  <a><img src="https://badges.aleen42.com/src/cli.svg" width="6.5%"/></a>
  <a><img src="https://img.shields.io/badge/Unix%20based%20os-053766?style=for-the-badge&logo=Linux&logoColor=white&" width="13%"/></a>
  <a><img src="https://img.shields.io/badge/Python3.7+-FFD43B?style=for-the-badge&logo=python&logoColor=blue" width="11.5%"/></a>
  <a><img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" width="4.6%"/></a>
  <a><img src="https://img.shields.io/badge/Docker-2CA5E0?style=for-the-badge&logo=docker&logoColor=white" width="8.5%"/></a>


  <a href="https://github.com/i2bc/SURFMAP/releases/latest" target="_blank"><img src="https://img.shields.io/github/release/i2bc/SURFMAP.svg?style=flat" alt="SURFMAP Release"/></a>
  <a href="https://opensource.org/licenses/MIT" target="_blank"><img src="https://img.shields.io/badge/License-MIT-yellow.svg?style=flat" alt="MIT License"/></a>
  <!-- <br> -->

  <!-- <img src="https://badgen.net/docker/pulls/lopesi2bc/surfmap?icon=docker&label=pulls" alt="Docker Pulls"/>
  <img src="https://img.shields.io/github/stars/i2bc/SURFMAP" alt="GitHub Stars"/>
  <img src="https://img.shields.io/github/forks/i2bc/SURFMAP.svg?style=flat" alt="GitHub Fork"/>
  <img src="https://img.shields.io/github/watchers/i2bc/SURFMAP.svg?style=social&style=plastic" alt="GitHub Watchers"/> -->
</div>

---



# SURFMAP

<div align="center">
  <img src="./doc/images/toc_Schweke_SURFMAP_cmyk.png" width="80%"/>  
</div>

# Table of contents

- [Aims](#Aims)
- [Installation](#Installation)
- [How it works](#How-it-works)
- [Usage of SURFMAP](#Usage-of-SURFMAP)
- [Supporting the project](#Supporting-the-project)
- [Contacts](#Contacts)
- [Licence](#Licence)
- [How to cite SURFMAP](#How-to-cite-SURFMAP)
- [References](#References)

# Aims
[Go to the top](#Table-of-contents)

<div>
<img src="./doc/images/TOC_Schweke_manuscript_revisions_forGitHub.png" width="60%" align="right"/>

SURFMAP is a free standalone and easy-to-use command-line interface (CLI) software that enables the fast and automated 2D projection of either predefined features of protein surface (electrostatic potential, Kyte-Doolittle hydrophobicity, Wimley-White hydrophobicity, stickiness and surface relief) or any descriptor encoded in the temperature factor column of a PDB file. The 2D maps computed by SURFMAP can be used to analyze and/or compare protein surface properties.
</div>


# Installation
[Go to the top](#Table-of-contents)

### Requirements

SURFMAP is a CLI tool that requires a UNIX-based OS system. It is written in python (version 3.7), R (version 3.6). It relies on the already included MSMS software ([1](#ref-1)) and may optionally require APBS ([2](#ref-2)) if the user wants to perform electrostatics calculations.

All those requirements (including APBS) are met in a [predefined Docker image](https://hub.docker.com/r/lopesi2bc/surfmap/tags) that we recommend the user to use. 

<details open>
<summary><b>For a usage of the docker image</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [Docker](https://docs.docker.com/get-docker/)

</details>

<details>
<summary><b>For a usage on your local OS</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [R >= 3.6](https://cran.r-project.org/)
- [APBS](https://github.com/Electrostatics/apbs/releases) (optional - only if you want to compute electrostatic potential)
 
</details>
<br>

> :bell: Please note that **whether you want to use the Docker image of SURFMAP or not, you will still need to [install the SURFMAP package](#How-to-install-SURFMAP)**. Indeed the package contains internal features that make the use of the Docker image totally transparent for the user who will not have to enter 'complex' commands for the connection of useful mounting points. In fact, the SURFMAP commands are almost exactly the same between the use of the docker image or not (see [here](#cmd_docker_or_not)).


### Recommendation

We strongly recommend that you install the SURFMAP package and its python dependencies in an isolated environment. Click in the section below for a short illustration on why and how to use an isolated environment.

<details style="margin-left: 32px">
<summary>How to use an isolated environment (recommended)</summary>
<br>
<p>
By using an isolated environment you will avoid potential version conflicts between python libraries when working on different projects. Some of the most popular tools to work with isolated python environments are [virtualenv](https://pypi.org/project/virtualenv/), [pyenv](https://pypi.org/project/pyenv/), [pipenv](https://pypi.org/project/pipenv/). 
</p>

Below is an example on how to use [virtualenv](https://pypi.org/project/virtualenv/).

#### 1. Install virtualenv
```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install virtualenv
python3 -m pip install virtualenv
```

#### 2. Create and activate an isolated environment
```bash
# create an isolated environment named 'myenv' (to adapt)
virtualenv myenv

# activate your isolated environment
source myenv/bin/activate
```

Once activated, any python library you'll install using pip will be installed in this isolated environment, and python will only have access to these packages.

Once you're done working on your project, simply type `deactivate` to exit the environment.
</details>


## How to install SURFMAP
[Go to the top](#Table-of-contents)

First, make sure you meet the [system requirements](#requirements) outlined earlier and consider the [recommendation](#recommendation). Then, follow instructions described in option 1 or 2 if you're not interested in accessing/modifying the source code, otherwise prefer option 3. 

<a id="install_option1"></a>
<details open>
<summary><h4>Option 1: from the archive (git not required)</h4></summary>

First download an archive of our latest release <a href="https://github.com/i2bc/SURFMAP/releases/latest" target="_blank">here</a>.

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP vx.x.x
python3 -m pip install SURFMAP-x.x.x.zip # (or .tar.gz) 
```
</details>


<details>
<summary><h4>Option 2: from the version control systems</h4></summary>

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP vx.x.x
python -m pip install -e git+https://github.com/i2bc/SURFMAP.git@v2.1.0#egg=surfmap
```
</details>

<details>
<summary><h4>Option 3: from this project repository</h4></summary>

```bash
# clone SURFMAP on your machine
git clone https://github.com/i2bc/SURFMAP.git

# go in the SURFMAP/ directory
cd SURFMAP

# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP
python3 -m pip install -e .
```
</details>


# How it works
[Go to the top](#Table-of-contents)


### SURFMAP workflow: inputs/outputs

<div align="center">
  <img src="./doc/images/surfmap_workflow.png" width="70%"/>

<i>The figure above represents the main steps of the SURFMAP worflow to compute the projection on a 2D map of a protein surface feature. More details about each step can be found in our article: see the [published version](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01269) or its [free version](https://www.biorxiv.org/content/10.1101/2021.10.15.464543v1)</i>
</div>
<br>

SURFMAP accepts as input either a *PDB file* or a *text file in a SURFMAP-specific matrix format*.
<br>
<br>

[Using a PDB file as input](#from-a-pdb-structure) is the most classic usage of SURFMAP. In this case, two outputs are generated: 
- the 2D map projection in a PDF format (PNG is also available)
- a matrix text file written in a SURFMAP-specific format

The matrix text file contains all information about each projected surface residue and their associated feature value. As the above figure shows, this text file is the direct input for the last step of the SURFMAP workflow as it is read to generate the 2D map projection.
<br>
<br>

[Using a text file in a SURFMAP-specific matrix format as input](#from-a-surfmap-matrix-file) represents a special case that could be useful if the user wants to generate a 2D map from an internally pre-processed matrix, such as to normalize or average with other matrices.

<details>
<summary>Example of a SURFMAP-specific matrix format (.txt)</summary>

<pre>absc    ord     svalue  residues
5       5       Inf     NA
5       10      Inf     NA
5       15      Inf     NA
...
5       80      Inf     GLU_120_A
5       85      Inf     GLU_120_A, GLN_301_A
5       90      Inf     GLN_301_A
5       95      Inf     GLN_301_A
5       100     Inf     GLN_301_A
5       105     Inf     GLN_301_A
...
360	175	Inf	NA
360	180	Inf	NA
</pre>
</details>


### Calling SURFMAP with Docker or not

Whether you want to use SURFMAP through a Docker or not, the commands are almost exactly the same. Indeed, in order to use the Docker image of SURFMAP, you will just have to add the CLI option **`--docker`**. If you want to use SURFMAP from an installation on your local OS, then simply remove this option. For example:

<a id="cmd_docker_or_not"></a>

```bash
# a command that will run on a Docker container
surfmap -pdb foo.pdb -tomap stickiness --docker

# the same command that will run on your local OS
surfmap -pdb foo.pdb -tomap stickiness
```

If the Docker image of SURFMAP is missing from your system, it will be automatically downloaded the first time you will execute a SURFMAP command.

> :bell: The version of the SURFMAP Docker image used is the same as the version of SURFMAP you will have installed. You can check your current version with the command `surfmap -v`. Yet if you want to use [another version of the SURFMAP Docker image](https://hub.docker.com/r/lopesi2bc/surfmap/tags), you will have to set a `SURFMAP_DOCKER_VERSION` environment variable with a value corresponding to an available tag version (e.g. `export SURFMAP_DOCKER_VERSION=2.1.0`).


# Usage of SURFMAP
[Go to the top](#Table-of-contents)

Once you have [installed the SURFMAP package](#how-to-install-surfmap), you should be ready to use SURFMAP. 

#### The example directory
To guide the user in the usage of SURFMAP, we will make use of files that you can find in the `example/` directory of SURFMAP. You can see where this directory is located on your machine with the following command:

```bash
python3 -c "import surfmap; print(surfmap.PATH_TO_EXAMPLES)"
```

Please note that for all command examples illustrated below, we will make [use of the Docker image of SURFMAP](#use-surfmap-with-docker-or-not).


#### SURFMAP options

<details>
<summary>List of all SURFMAP options</summary>

<pre>usage: surfmap [-h] (-pdb PDB | -mat MAT | -v) -tomap TOMAP [-proj PROJ] [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth] [--png] [--keep]
               [--docker] [--pqr PQR] [-ff FF] [-verbose VERBOSE]

options:
  -h, --help        show this help message and exit
  -pdb PDB          Path to a PDB file
  -mat MAT          Input matrix. If the user gives an imput matrix, SURFMAP will directly compute a map from it.
  -v, --version     Print the current version of SURFMAP.
  -tomap TOMAP      Specific key of the feature to map. One of the following: stickiness, kyte_doolittle, wimley_white, electrostatics,
                    circular_variance, bfactor, binding_sites, all.
  -proj PROJ        Choice of the projection. Argument must be one of the following: flamsteed, mollweide, lambert. Defaults to flamsteed.
  -res RES          File containing a list of residues to map on the projection. Expected format has the following space/tab separated column
                    values: chainid resid resname
  -rad RAD          Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the
                    smoother the surface. Defaults to 3.0
  -d D              Output directory where all files will be written. Defaults to &apos;./output_SURFMAP_$pdb_$tomap&apos; with $pdb and $tomap based on
                    -pdb and -tomap given values
  -s S              Value defining the size of a grid cell. The value must be a multiple of 180. Defaults to 5.0.
  --elec-max-value ELEC_MAX_VALUE
                        Maximum value to be used for the electrostatics color scale. The value will be converted as an absolute value to make the scale symetric around 0. For
                        instance, a value of 5.63 will scale the electrosctatics color values from -5.63 to 5.63.
  --bfactor-min-value BFACTOR_MIN_VALUE
                        Minimum value to be used for the bfactor color scale.
  --bfactor-max-value BFACTOR_MAX_VALUE
                        Maximum value to be used for the bfactor color scale.
  --nosmooth        If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)
  --png             If chosen, a map in png format is computed (default: only pdf format is generated)
  --keep            If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)
  --docker          If chosen, SURFMAP will be run on a docker container (requires docker installed).
  --pqr PQR         Path to a PQR file used for electrostatics calculation. Option only available if &apos;-tomap electrosatics&apos; is requested.
                    Defaults to None.
  -ff FF            Force-field used by pdb2pqr for electrostatics calculation. One of the following: AMBER, CHARMM, PARSE, TYL06, PEOEPB,
                    SWANSON. Defaults to CHARMM.
  -verbose VERBOSE  Verbose level of the console log. 0 for silence, 1 for info level, 2 for debug level. Defaults to 1.
</pre>
</details>


## Projection of a protein surface feature on a 2D map

In order to generate a 2D map projection of a protein surface feature, two inputs are required:
- either a PDB file (`-pdb` option) OR a matrix text file written in a SURFMAP-specific format (`-mat` option)
- a valid key referring to a feature to map (listed in the table below)

| Valid feature key | Feature details |
| --- | --- |
`kyte_doolittle` | Residue hydrophobicity directly derived from the Kyte-Doolittle scale ([3](#ref-3))
`wimley_white` | Residue hydrophobicity directly derived from the Wimley-White scale ([4](#ref-4))
`stickiness` | Propensity of each amino acid to be involved in protein−protein interfaces ([5](#ref-5))
`circular_variance` | Descriptor of the local (residue scale) geometry of a surface region: low values reﬂects protruding residues, while high values indicates residues located in cavities ([6](#ref-6))
`circular_variance_atom` | Descriptor of the local geometry (atomic scale) of a surface region: low values reﬂects protruding atoms, while high values indicates atoms located in cavities. ([6](#ref-6))
`electrostatics` | Electrostatic potential of the protein surface (atomic scale) - Requires the APBS software ([2](#ref-2))
`bfactor` | Any feature stored in the temperature factor of the input PDB ﬁle
`all` | Compute sequentially the following features: `kyte_doolittle`, `wimley_white`, `stickiness` and `circular_variance`

#### From a PDB structure

```bash
# example - command to map the stickiness values for residues at the surface of the chain A of 1g3n.pdb
surfmap -pdb 1g3n_A.pdb -tomap stickiness --docker
```

The output has the following structure and content:
<pre><font color="#12488B"><b>output_SURFMAP_1g3n_A_stickiness/</b></font>
├── <font color="#12488B"><b>maps</b></font>
│   └── 1g3n_A_stickiness_map.pdf
├── parameters.log
├── surfmap.log
└── <font color="#12488B"><b>smoothed_matrices</b></font>
    └── 1g3n_A_stickiness_smoothed_matrix.txt
</pre>

with:
- `parameters.log`: a summary of the parameters used to compute the map
- `surfmap.log`: a log file of each of the step of the SURFMAP workflow
- `1g3n_A_stickiness_map.pdf`: the generated 2D map in PDF format
- `1g3n_A_stickiness_smoothed_matrix.txt`: a computed smoothed matrix file (txt file) used to generate the 2D map. This matrix has the expected format of a matrix file that can be used as a direct input of SURFMAP through the used of the `-mat` argument.


<details>
<summary><h5>Note on electrostatics calculations</h5></summary>

The electrostatics potential is calculated through the use of APBS and is initially based on the generation of a PQR file which will contain the charge and radius of each atom in the input PDB file. In SURFMAP, this PQR file is generated through the use of [pdb2pqr](https://pdb2pqr.readthedocs.io/en/latest/) which reads atomic parameters from a force field accessible from its package itself. While the CHARMM force field is used by default in SURFMAP, all force fields accessible in [pdb2pqr](https://pdb2pqr.readthedocs.io/en/latest/) (AMBER, CHARMM, PARSE, TYL06, PEOEPB, SWANSON) can be used in SURFMAP with the `-ff` option. For example:

```bash
# will use the CHARMM force-field (default)
surfmap -pdb 1g3n_A.pdb -tomap electrostatics --docker

# will use the AMBER force-field
surfmap -pdb 1g3n_A.pdb -tomap electrostatics -ff AMBER --docker
```

For the particular case where a user would like to compute electrostatics potential with any other force-field (e.g. for a coarse-grained PDB file), SURFMAP can be used with the additional option `-pqr` that must be followed with a PQR file generated by the user himself. For example:

```bash
# will read atomic parameters from the PQR given as input
surfmap -pdb structure.pdb -tomap electrostatics -pqr structure.pqr --docker
```

</details>


#### From a SURFMAP matrix file

A matrix written in a SURFMAP-specific format can also be used as an input to generate a 2D map. The feature to map has to be the same as the one used to generate the matrix file. As a fancy usage example, the command below will reproduce the 2D map generated from the command above:

```bash
# example - command to create a map from a SURFMAP matrix file generated with stickiness values
surfmap -mat output_SURFMAP_1g3n_A_stickiness/smoothed_matrices/1g3n_A_stickiness_smoothed_matrix.txt -tomap stickiness --docker
```

A more realistic usage of this option would be to compute maps from your internally pre-processed matrices. For example you may have generated 2D maps of a same protein in different conformational states and then may want to compute an averaged matrix file (please note that we don't provide such script utilities).

<details>
<summary>Example of a SURFMAP-specific matrix format (.txt)</summary>

<pre>absc    ord     svalue  residues
5       5       Inf     NA
5       10      Inf     NA
5       15      Inf     NA
...
5       80      Inf     GLU_120_A
5       85      Inf     GLU_120_A, GLN_301_A
5       90      Inf     GLN_301_A
5       95      Inf     GLN_301_A
5       100     Inf     GLN_301_A
5       105     Inf     GLN_301_A
...
360	175	Inf	NA
360	180	Inf	NA
</pre>
</details>



## Projection of interface residues on a 2D map

<p><img src="./doc/custom/project_binding_site.svg"></p>

Instead of projecting a protein surface feature on a 2D map, you may be interested in the projection of interface residues. This is possible with the option `-tomap binding_sites` of SURFMAP. 

With the `-tomap binding_sites` option, a discrete color scale is used to associate one color to each different value found in the b-factor column. So in order to use this option, your input PDB file must contain discrete values in the b-factor column for each atoms, the value depending on whether the atoms belong to an interface or not. For example:
- `0` for atoms that are not part of any binding sites
- `1` for atoms being part of one known binding site
- `2` for atoms being part of a second binding site (if there is)
- `...`


We provide two utility scripts to help users generating a PDB file that can be used with the `-tomap binding_sites` option of SURFMAP:
- `extract_interface`
- `write_pdb_bs`

##### Usage of extract_interface

From multi-chain PDB file, the command `extract_interface` will find the interface residues between a given chain (or set of chains) and all the other chains of the input PDB structure. It will then output a new PDB file of the given chain(s) with the expected format for the `-tomap binding_sites` option.

The command below illustrates the usage of `extract_interface` with the PDB file `1g3n_ABC.pdb` in the example directory. 

```bash
# generate a PDB file of the chain A in which the b-factor column will contain a discrete value for each different interface residues that will be found between chains A and B, and chains A and C
extract_interface -pdb 1g3n_ABC.pdb -chains A
```

It will generate two output files:
- `1g3n_ABC_chain-A_bs.pdb`: a PDB file ready for use by the command `surfmap` with the option `-tomap binding_sites`.
- `1g3n_ABC_chain-A_interface.txt`: a text file containing information about identified interface residues. This file can be edited and used as input for the command `write_pdb_bs` described below.


So now, we can map interface residues of the chain A of 1G3N:
```bash
# Use the PDB file generated with the command above to project labelled residues on a 2D map 
surfmap -pdb 1g3n_ABC_chain-A_bs.pdb -tomap binding_sites --docker
```

##### Usage of write_pdb_bs

The command `write_pdb_bs` is made to avoid the manual editing of the b-factor column of a PDB file that you would like to use with the `-tomap binding_sites` option. The command takes as inputs:
- a PDB file of interest
- a text file listing interface residues of interest

The text file listing interface residues must be formatted as follows:
- 1st column: interface residue chain name
- 2nd column: residue ID
- 3rd column: residue name
- 4th column: discrete value (one value per different binding site; one color will be assigned to each different value)

<details>
<summary>Example of a text file listing interface residues</summary>

<pre>A	14	GLU	1
A	15	CYS	1
A	16	VAL	1
...
A	155	SER	1
A	156	SER	1
A	47	VAL	2
A	49	THR	2
A	50	GLY	2
...
A	139	HIS	2
A	140	ARG	2
A	292	TYR	2

</pre>
</details>


As a fancy example, the command below will reproduce the PDB file `1g3n_ABC_chain-A_bs.pdb` ready for use by `surfmap` with the option `-tomap binding_sites`:
```bash
write_pdb_bs -pdb 1g3n_ABC_chain-A_bs.pdb -res 1g3n_ABC_chain-A_interface.txt
```

The output file will have the basename of the PDB file given as input with the suffix `_bs.pdb`


# Supporting the project
[Go to the top](#Table-of-contents)

- If you find a bug or have a suggestion for a new feature, please report it via an [issue](https://github.com/i2bc/SURFMAP/issues)
- If you find SURFMAP useful, consider starring the repository


# Contacts
[Go to the top](#Table-of-contents)

If you have any question regarding SURFMAP, you can contact us:
- [@anne.lopes@i2bc.paris-saclay.fr](mailto:anne.lopes@i2bc.paris-saclay.fr) (project leader)
- [@hugo.schweke@weizmann.ac.il](mailto:hugo.schweke@weizmann.ac.il) (original code author)
- [@nicolas.chevrollier@laposte.net](mailto:nicolas.chevrollier@laposte.net) (project maintainer)


# Licence
[Go to the top](#Table-of-contents)

This project is under the MIT License terms. Please have a look at the LICENSE file for more details.


# How to cite SURFMAP
[Go to the top](#Table-of-contents)

If SURFMAP has been useful to your research, please cite us as well as the original MSMS paper:

> Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier, Simon Gosset, Anne Lopes. SURFMAP: a software for mapping in two dimensions protein surface features. J. Chem. Inf. Model. 2022. [Link](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01269)

> Sanner, M. F., Olson A.J. & Spehner, J.-C. (1996). Reduced Surface: An Efficient Way to Compute Molecular Surfaces. Biopolymers 38:305-320. [Link](https://doi.org/10.1002/(sici)1097-0282(199603)38:3%3c305::aid-bip4%3e3.0.co;2-y)
<br>

Moreover, if you use APBS in your research, please cite one or more of the following papers listed in the [Supporting APBS](https://apbs.readthedocs.io/en/latest/supporting.html) documentation page.
<br>

# References
[Go to the top](#Table-of-contents)

<a id="ref-1"></a>

> (1) Michel Sanner, Arthur J. Olson, Jean Claude Spehner (1996). Reduced Surface: an Efficient Way to Compute Molecular Surfaces. Biopolymers, Vol 38, (3), 305-320.


<a id="ref-2"></a>

> (2) Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112-128, 2018.


<a id="ref-3"></a>

> (3) Kyte, J.; Doolittle, R. F. A Simple Method for Displaying the Hydropathic Character of a Protein. J. Mol. Biol. 1982, 157 (1), 105−132.


<a id="ref-4"></a>

> (4) Wimley, W. C.; White, S. H. Experimentally Determined Hydrophobicity Scale for Proteins at Membrane Interfaces. Nat. Struct. Biol. 1996, 3 (10), 842−848.


<a id="ref-5"></a>

> (5) Levy, E. D. A Simple Definition of Structural Regions in Proteins and Its Use in Analyzing Interface Evolution. J. Mol. Biol. 2010, 403 (4), 660−670.


<a id="ref-6"></a>

> (6) Mezei, M. A New Method for Mapping Macromolecular Topography. J. Mol. Graph. Model 2003, 21 (5), 463−472.
