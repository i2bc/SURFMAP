<div align="center">

  <a><img src="https://badges.aleen42.com/src/cli.svg" width="6.5%"/></a>
  <a><img src="https://img.shields.io/badge/Unix%20based%20os-053766?style=for-the-badge&logo=Linux&logoColor=white&" width="13%"/></a>
  <a><img src="https://img.shields.io/badge/Python3.7+-FFD43B?style=for-the-badge&logo=python&logoColor=blue" width="11.5%"/></a>
  <a><img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" width="4.6%"/></a>
  <a><img src="https://img.shields.io/badge/Docker-2CA5E0?style=for-the-badge&logo=docker&logoColor=white" width="8.5%"/></a>


  <img src="https://img.shields.io/github/release/i2bc/SURFMAP.svg?style=flat" alt="GitHub Release"/>
  <img src="https://img.shields.io/github/license/i2bc/SURFMAP.svg?style=flat" alt="GitHub License"/>
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
- [How it works](#How-it-works)
- [Installation](#Installation)
- [Usage of SURFMAP](#Usage-of-SURFMAP)
- [How to cite SURFMAP](#How-to-cite-SURFMAP)

# Aims
[Go to the top](#Table-of-contents)

<div>
<img src="./doc/images/TOC_Schweke_manuscript_revisions_forGitHub.png" width="60%" align="right"/>

SURFMAP is a free standalone and easy-to-use command-line interface (CLI) software that enables the fast and automated 2D projection of either predefined features of protein surface (electrostatic potential, Kyte-Doolittle hydrophobicity, Wimley-White hydrophobicity, stickiness and surface relief) or any descriptor encoded in the temperature factor column of a PDB file. The 2D maps computed by SURFMAP can be used to analyze and/or compare protein surface properties.
</div>


# How it works
[Go to the top](#Table-of-contents)

<div align="center">
  <img src="./doc/images/surfmap_workflow.png" width="75%"/>

<i>The figure above represents the main steps of the SURFMAP worflow to compute the projection on a 2D map of a protein surface feature. More details about each step can be found in [our article](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01269).</i>
</div>

<br>


SURFMAP accepts as input either a *PDB file* or a *text file in a SURFMAP-specific matrix format*

[Using a PDB file as input](#from-a-pdb-structure--pdb) is the most classic usage of SURFMAP. In this case, two outputs are generated: the 2D map projection in a PDF format (PNG is also available) and a text file. This text file is a matrix written in a SURFMAP-specific format containing all information about each projected surface residue and their associated feature value. As the above figure shows, this text file is the direct input for the last step of the SURFMAP workflow as it is read to generate the 2D map projection.

<details>
<summary><b>Example of a SURFMAP-specific matrix format (.txt)</b></summary>

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

[Using a text file in a SURFMAP-specific matrix format as input](#from-a-surfmap-matrix-file--mat) to SURFMAP represents a special case that could be useful if the user wants to generate a 2D map from an internally pre-processed matrix, such as to normalize or average with other matrices.


# Installation
[Go to the top](#Table-of-contents)

SURFMAP is a CLI tool that requires a UNIX-based OS system. It is written in python (version 3.7), R (version 3.6). It relies on the already included MSMS software (1) and may optionally require APBS (2) if the user wants to perform electrostatics calculations.

All those requirements (including APBS) are met in a [predefined Docker image](https://hub.docker.com/r/lopesi2bc/surfmap) that we recommend the user to use. 

**Please also note that whether you want to use the Docker image of SURFMAP or not, you will still need to [install the SURFMAP package](#install_option1)**. Indeed the package contains internal features that make the use of the Docker image totally transparent for the user who will not have to enter 'complex' commands for the connection of useful mounting points. In fact, the SURFMAP commands are almost exactly the same between the use of the docker image or not (see [here](#cmd_docker_or_not)).

See below the requirements that must be met for the use of SURFMAP through its Docker image or a complete local installation.

### Requirements

<details>
<summary><b>For a usage of the docker image</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [Docker](https://docs.docker.com/get-docker/)
- the SURFMAP docker image: `docker pull lopesi2bc/surfmap`

</details>

<details>
<summary><b>For a complete local installation</b></summary>

- an UNIX-based OS system (any linux distribution, a MacOS system or [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) on windows)
- [Python >= 3.7](https://www.python.org/downloads)
- [R >= 3.6](https://cran.r-project.org/)
- [APBS](https://github.com/Electrostatics/apbs/releases) (optional - only if you want to compute electrostatic potential)
 
</details>

#### Notes

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
Choose option 1 or 2 if you're not interested in accessing/modifying the source code, otherwise prefer option 3. 

<a id="install_option1"></a>
<details open>
<summary><h4>Option 1: from the archive (git not required)</h4></summary>

First download an archive of our latest release <a href="https://github.com/i2bc/SURFMAP/releases/latest" target="_blank">here</a>.

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP v2.0.0
python3 -m pip install SURFMAP-v2.0.0.zip # (or .tar.gz) 
```
</details>


<details>
<summary><h4>Option 2: from the version control systems</h4></summary>

```bash
# upgrade pip to its latest version
python3 -m pip install --upgrade pip

# install SURFMAP v2.0.0
python -m pip install -e git+https://github.com/i2bc/SURFMAP.git@v2.0.0#egg=surfmap
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


# Usage of SURFMAP
[Go to the top](#Table-of-contents)


#### The example directory
To guide the user in the usage of SURFMAP, we will make use of files that you can find in the `example/` directory of SURFMAP. You can see where this directory is located on your machine with the following command:

```bash
python3 -c "import surfmap; print(surfmap.PATH_TO_EXAMPLES)"
```

All command examples will make use of the docker image of SURFMAP thanks to the CLI option **`--docker`**. If you want to use SURFMAP through a local install, then simply remove this option. For instance:

<a id="cmd_docker_or_not"></a>

```bash
# a command that will run on a docker container
surfmap -pdb foo.pdb -tomap stickiness --docker

# the same command that will run locally
surfmap -pdb foo.pdb -tomap stickiness
```


## Projection of a protein surface feature on a 2D map

In order to generate a 2D map, two inputs are required:
- either a PDB file (`-pdb` option) OR a SURFMAP matrix file (`-mat` option)
- a valid key referring to a feature to map:
  - `kyte_doolittle`
  - `wimley_white`
  - `stickiness`
  - `circular_variance`
  - `circular_variance_atom`
  - `electrostatics` (requires APBS)
  - `bfactor`
  - `all` (to compute sequentially kyte_doolittle, wimley_white, stickiness and circular_variance features)

### From a PDB structure `-pdb`

```bash
# example - command to map the stickiness values for residues at the surface of the chain A of 1g3n.pdb
surfmap -pdb 1g3n_A.pdb -tomap stickiness --docker
```

The output has the following structure and content:
<pre><font color="#12488B"><b>output_SURFMAP_1g3n_A_stickiness/</b></font>
├── <font color="#12488B"><b>maps</b></font>
│   └── 1g3n_A_stickiness_map.pdf
├── parameters.log
└── <font color="#12488B"><b>smoothed_matrices</b></font>
    └── 1g3n_A_stickiness_smoothed_matrix.txt
</pre>

with:
- `parameters.log`: a summary of the parameters used to compute the map
- `1g3n_A_stickiness_map.pdf`: the generated 2D map in PDF format
- `1g3n_A_stickiness_smoothed_matrix.txt`: a computed smoothed matrix file (txt file) used to generate the 2D map. This matrix has the expected format of a matrix file that can be used as a direct input of SURFMAP through the used of the `-mat` argument.


### From a SURFMAP matrix file `-mat`

A matrix written in a SURFMAP-specific format can also be used as an input to generate a 2D map. The feature to map has to be the same as the one used to generate the matrix file. As a fancy usage example, the command below will reproduce the 2D map generated from the command above:

```bash
# example - command to create a map from a SURFMAP matrix file generated with stickiness values
surfmap -mat output_SURFMAP_1g3n_A_stickiness/smoothed_matrices/1g3n_A_stickiness_smoothed_matrix.txt -tomap stickiness --docker
```

A more realistic usage of this option is to compute maps from your own customized matrices. For example you can create maps of a same protein in different conformational states. In this case, you may want to compute an averaged matrix file (please note that we don't provide such script utilities).


## Projection of interface residues on a 2D map

<p><img src="./doc/custom/project_binding_site.svg"></p>

Instead of projecting a protein surface feature on a 2D map, you may be interested in the projection of interface residues. This is possible with the option `-tomap binding_sites` of SURFMAP. 

With the `-tomap binding_sites` option, a discrete color scale is used to associate one color to each different value found in the b-factor column. So in order to use this option, your input PDB file must contain discrete values in the b-factor column for each atoms, the value depending on whether the atoms belong to an interface or not. For instance:
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
# generate a PDB file of the chain A in which the b-factor column will contain a discrete value for each different interface residues that will be found between chains A and B, and A and C
extract_interface -pdb 1g3n_ABC.pdb -chains A
```

It will generate two output files:
- `1g3n_ABC_chain-A_bs.pdb`: a PDB file ready for use by the command `surfmap` with the option `-tomap binding_sites`.
- `1g3n_ABC_chain-A_interface.txt`: a text file containing information about identified interface residues. This file can be edited and used as input for the command `write_pdb_bs` described below.


So now, we can map interface residues of the chain A of 1G3N:
```bash
# Use the PDB file generated with the command above to project labelled residues on a 2D map 
surfmap -pdb 1g3n_ABC_chain-A_bs.pdb -tomap binding_sites
```

##### Usage of write_pdb_bs

The command `write_pdb_bs` is made to avoid the manual editing of the b-factor column of a PDB file that you would like to use with the `-tomap binding_sites` option. The command takes as inputs:
- a PDB file of interest
- a text file listing interface residues of interest

The text file listing interface residues must be formatted as follows:
- 1st column: the chain name that has an interface residue
- 2nd column: a residue ID
- 3rd column: a residue name
- 4th column: a discrete value (one value per different binding site)

For instance, if you want to map the residues GLU-14 and CYS-15 of the chain A as part of an interface, and the residues GLY-50 and GLU-51 of the chain A as part of another interface, you should create your input file as follows:
```
A	14	GLU	1
A	15	CYS	1
A	50	GLY	2
A	17	GLU	2
```

As a fancy example, the command below will reproduce the PDB file `1g3n_ABC_chain-A_bs.pdb` ready for use by `surfmap` with the option `-tomap binding_sites`:
```bash
write_pdb_bs -pdb 1g3n_ABC_chain-A_bs.pdb -res 1g3n_ABC_chain-A_interface.txt
```

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

(1) Michel Sanner, Arthur J. Olson, Jean Claude Spehner (1996). Reduced Surface: an Efficient Way to Compute Molecular Surfaces. Biopolymers, Vol 38, (3), 305-320.

(2) Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112-128, 2018.
