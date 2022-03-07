# SURFMAP
<p align="center">
  <img src="./doc/images/toc_Schweke_SURFMAP_cmyk.png" width="80%"/>
</p>

# Table of contents

- [Aims](#Aims)
- [Prerequisites](#Prerequisites)
- [Download](#Download)
- [Usage of SURFMAP](#Usage-of-SURFMAP)
- [Use of the pre-built docker image of SURFMAP](#Use-of-the-pre-built-docker-image-of-SURFMAP)
- [Installing APBS](#Installing-APBS)
- [How to cite SURFMAP](#How-to-cite-SURFMAP)


# Aims
[Go to the top](#Table-of-contents)

<div>
<img src="./doc/images/toc_Schweke_SURFMAP_rev.png" width="60%" align="right"/>

SURFMAP is a free standalone and easy-to-use software that enables the fast and automated 2-D projection of either predefined features of protein surface (electrostatic potential, Kyte-Doolittle hydrophobicity, Wimley-White hydrophobicity, stickiness and surface relief) or any descriptor encoded in the temperature factor column of a PDB file. The 2-D maps computed by SURFMAP can be used to analyze and/or compare protein surface properties.
</div>


# Prerequisites
[Go to the top](#Table-of-contents)

SURFMAP is written in python (version 3.7), R (version 3.6) and bash. It relies on the MSMS software (1) and may optionally requires APBS (2) if the user wants to perform electrostatics calculations.
All these requirements (including the APBS software) are fullfilled in a [pre-built docker image of SURFMAP](#Use-of-the-pre-built-docker-image-of-SURFMAP).

**Please take note that we strongly advise to use the docker image, especially if you want to compute electrostatics potential. Indeed the installation process of APBS can be tricky. By using the pre-built docker image of SURFMAP, you will not have to install anything except docker itself.**

If you don't want to use the pre-built docker image of SURFMAP, you will need:
- an UNIX-based OS system (any linux distribution, a MacOS system or WSL2 on windows)
- python >= 3.7
- R >= 3.6

# Download
[Go to the top](#Table-of-contents)

The recommended way to retrieve the project is to clone the repository.

To specifically clone this branch type in a terminal:
```
git clone -b feature/new_install https://github.com/i2bc/SURFMAP.git
```

Otherwise type the following to clone the master branch (by default):
 with the following command:
```
git clone https://github.com/i2bc/SURFMAP.git
```
It will allow you to easily update SURFMAP with its latest version through the command `git pull`




Alternatively, you can click <a href="https://github.com/i2bc/SURFMAP/releases/latest" target="_blank">here</a> to access the latest source code of the release.

After downloading, you'll need to unarchive the project. If you downloaded:
- the .zip file: type on a terminal (Unix and macOS): `unzip SURFMAP-x.x.zip`
- the .tar.gz file: type on a terminal (Unix and macOS): `tar xzvf SURFMAP-x.x.tar.gz`

Once downloaded, please go to the SURFMAP directory (`cd SURFMAP/` where you should have the following files/directories:

<pre><font color="#3465A4"><b>.</b></font>
├── <font color="#3465A4"><b>doc</b></font>
├── <font color="#3465A4"><b>example</b></font>
├── how_to_cite_us.txt
├── LICENSE
├── <font color="#3465A4"><b>MSMS</b></font>
├── README.md
├── requirements.txt
├── <font color="#4E9A06"><b>run_surfmap.py</b></font>
├── <font color="#3465A4"><b>scripts</b></font>
├── <font color="#4E9A06"><b>SURFMAP_launcher.py</b></font>
├── <font color="#3465A4"><b>tools</b></font>
└── <font color="#3465A4"><b>visualizer</b></font>
</pre>

## Install required python libraries
Simply type in a terminal:
```bash
pip3 install -r requirements.txt
```

# Usage of SURFMAP
[Go to the top](#Table-of-contents)

If you followed the steps above, SURFMAP should be ready to use through the script `SURFMAP_launcher.py`. To make it sure, type in a terminal the following command:
```
python3 SURFMAP_launcher.py
```
This command should display:

<pre>SURFMAP:    Projection of protein surface features on 2D map
Authors:    Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier,
            Simon Gosset, Anne Lopes
Version:    x.x (latest)
Copyright (c) 2021, H. Schweke

...
    
usage: SURFMAP_launcher.py [-h] -pdb PDB -tomap
                           {bfactor,kyte_doolittle,binding_sites,all,circular_variance_atom,stickiness,wimley_white,electrostatics,circular_variance}
                           [-proj {sin,aitoff,moll,cyl}] [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth] [--png]
                           [--keep]
SURFMAP_launcher.py: error: the following arguments are required: -pdb, -tomap
</pre>

To guide the user on how to use SURFMAP, we will use files in the `example/` directory that can be found in the downloaded SURFMAP project:

<pre><font color="#3465A4"><b>example/</b></font>
├── 1g3n_A.pdb
├── 1gv3_A_binding_sites.pdb
├── <font color="#3465A4"><b>example_1g3n_A_stickiness/</b></font>
├── <font color="#3465A4"><b>example_1g3n_A_stickiness_allfiles/</b></font>
├── <font color="#3465A4"><b>example_1g3n_A_stickiness_mapping_residues/</b></font>
├── <font color="#3465A4"><b>example_1gv3_A_stickiness_binding_sites/</b></font>
├── <font color="#3465A4"><b>example_mapping_stickiness/</b></font>
├── README
└── residues_to_map.txt
</pre>

## SURFMAP inputs and outputs

SURFMAP allows to compute different protein surface features and to map them on a 2-D plan through a sinusoidal projection by default (other projections are available: mollweide, aitoff and cylindric). Thus two mandatory arguments must be given as inputs by the user: `-pdb` and `-tomap`:

- the `-pdb` argument must be followed by the protein structure in PDB format the user wants to analyse
- the `-tomap` argument must be given a keyword representing the protein surface feature the user wants to map. The user can also use the option `all` to map the Kyte-Doolittle hydrophobicity, the Wimley-White hydrophobicity, the stickiness and the circular variance per residue at the same time. The available keywords are listed below (see SURFMAP_manual.pdf in `doc/` or the original article for a description):
  - wimley_white
  - kyte_doolittle
  - stickiness
  - circular_variance
  - circular_variance_atom
  - electrostatics
  - bfactor
  - binding_sites
  - all

For instance, the following command line will map the stickiness protein surface feature for the chain A of the protein [1G3N](https://www.rcsb.org/structure/1G3N):

```bash
python3 SURFMAP_launcher.py -pdb 1g3n_A.pdb -tomap stickiness
```

The above command will generate an output directory named `output_SURFMAP_1g3n_A_stickiness/` with the following content:

<pre><font color="#3465A4"><b>output_SURFMAP_1g3n_A_stickiness/</b></font>
├── log_parameters
├── <font color="#3465A4"><b>maps/</b></font>
│   └── 1g3n_A_stickiness_map.pdf
└── <font color="#3465A4"><b>smoothed_matrices/</b></font>
    └── 1g3n_A_stickiness_smoothed_matrix.txt
</pre>

with:
- `log_parameters`: a summary of the basic parameters used to compute the map
- `1g3n_A_stickiness_sinusoidal_map.pdf`: the generated 2-D map in pdf format
- `1g3n_A_stickiness_smoothed_matrix.txt`: a computed smoothed matrix file used to generate the 2-D map


## SURFMAP optional parameters

The following table lists the optional parameters that can be used when running SURFMAP:

| Optional parameters | Description |
| --- | --- |
| -proj | Choice of the projection. Argument must be one of the following: `sin` for sinusoidal, `moll` for mollweide, `aitoff` aitoff or `cyl` for cylindric equal area (default: `sin`) |
| -res | File containing a list of residues to map on the projection |
| -rad | Radius added to the usual atomic radius used to calculate the solvent excluded surface. The higher the radius the smoother the surface (default: 3.0 Angström) |
| -d | Output directory where all files will be written (default: './output_SURFMAP_$pdb_$tomap' where $pdb and $tomap are the inputs given to `-pdb` and `-tomap` arguments, respectiveley) |
| -s | Size of a grid cell. Necessary that 180%cellsize == 0 (default: 5) |
| --nosmooth | If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!) |
| --png | If chosen, a map in png format is computed (default: only pdf format is generated) |
| --keep | If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept) |

<br>

# Visualize interactively your maps

All maps generated with SURFMAP can be interactively inspected with the help of an `index.html` file located in `SURFMAP/visualizer/`.

To interactively visualize a map:
- open `index.html` in your favorite browser
- upload a smoothed matrix file (in the directory `smoothed_matrices` from an output generated with SURFMAP)

Once the map appears, you can click on any pixel to see the corresponding residue(s).



# Use of the pre-built docker image of SURFMAP
[Go to the top](#Table-of-contents)

This image will work on any system that can run docker (Unix, MacOS, or Windows 10 through WSL2). It includes all the dependencies and external softwares (MSMS & APBS) to run SURFMAP so that the user don't need to install anything on its machine (except docker).  

## 1. Get docker on your machine

You’ll first need to create an account on [docker hub](https://hub.docker.com/) and [install docker](https://docs.docker.com/get-docker/) on your machine.


## 2. Run SURFMAP with `run_surfmap.py`

Once you have docker installed on your machine, you are ready to use SURFMAP. Simply type in a terminal:
```bash
python3 run_surfmap.py
```  

The first time you’ll type it you may have to also type docker login in the terminal and fill the fields with your docker ID and docker password. It will then take few minutes to download the image on your machine. Once the SURFMAP image has been successfully downloaded on your machine, you should see after typing the `python3 run_surfmap.py` command in a terminal the following message:

<pre>
<font color="#4E9A06"><b>tutor@surfmap</b></font>:<font color="#3465A4"><b>~/i2bc/SURFMAP</b></font>$ python3 run_surfmap.py
usage: run_surfmap.py [-h] -pdb PDB -tomap
                      {electrostatics,all,circular_variance,circular_variance_atom,wimley_white,stickiness,kyte_doolittle,binding_sites,bfactor}
                      [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth]
                      [--png] [--keep]
run_surfmap.py: error: the following arguments are required: -pdb, -tomap
</pre>


# Tip

<a id="surfmap-alias"></a>
Wether you call SURFMAP through `SURFMAP_launcher.py` or `run_surfmap.py`, you can make those scripts callable from anywhere on your machine. To do so, add the following lines at the end of your `~/.bashrc` (or `~/.bash_profile` or `~/.profile`) file:

```bash
# Alias to run surfmap from its docker image"
alias surfmap='python3 PATH/TO/SURFMAP/run_surfmap.py'
```
or 
```bash
alias surfmap='python3 PATH/TO/SURFMAP/SURFMAP_launcher.py'
```

Make sure to replace `PATH/TO/SURFMAP/` with the absolute path you’ve downloaded SURFMAP to. Then type `source ~/.bashrc` (or `~/.bash_profile` or `~/.profile`) in the terminal to make the change effective.

From now on, you should be able to use SURFMAP by simply typing in a terminal `surfmap`.


# Installing APBS

**Please take note that we strongly advise to use the docker image, especially if you want to compute electrostatics potential. Indeed the installation process of APBS can be tricky. By using the pre-built docker image of SURFMAP, you will not have to install anything except docker itself.**


If you choose not to use the pre-built docker image of SURFMAP and want to compute electrostatics potential, you will also need to install APBS. We recommend to install the pre-compiled binaries from the version 3.0.0 that can be found [here](https://github.com/Electrostatics/apbs/releases/tag/v3.0.0).
The install documentation from the pre-compiled binaries is accessible [there](https://apbs.readthedocs.io/en/latest/getting/index.html#installing-from-pre-compiled-binaries).


After downloading the pre-compiled binaries, you will need to edit your `~/.bashrc` (or `~/.bash_profile` or `~/.profile`) file to make your system aware of the location of some APBS-3.0.0 required paths. Let's say you have downloaded the (Linux) pre-compiled binaries of APBS here: `~/APBS-3.0.0.Linux`. Then open your favorite text editor and add the following lines:

```bash
export PATH="$HOME/APBS-3.0.0.Linux/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/APBS-3.0.0.Linux/lib"
export APBS="$HOME/APBS-3.0.0.Linux/"
```

Finally type `source ~/.bashrc` (or `~/.bash_profile` or `~/.profile`) in the terminal. From now on the `apbs` command should be accessible.

You now have to fix a bug in a file used by APBS. To do so, simply type the following command:
```bash
sed -i '152s/nsmall/nsmall\[i\]/' $APBS/tools/manip/psize.py
```

Now SURFMAP should be ready for use. After typing in a terminal `apbs` you should see:

<pre>----------------------------------------------------------------------
    APBS -- Adaptive Poisson-Boltzmann Solver Version 3.0
    
    ...
    ...

    apbs [options] apbs.in

    where apbs.in is a formatted input file and [options] are:

--output-file=&lt;name&gt;     Enables output logging to the path
    listed in &lt;name&gt;.  Uses flat-file
    format is --output-format is not used.
--output-format=&lt;type&gt;   Specifies format for logging.  Options
    for type are either &quot;xml&quot; or &quot;flat&quot;.
--help                   Display this help information.
--version                Display the current APBS version.
----------------------------------------------------------------------
</pre>

### Notes on possible missing libraries

Depending on your system, you may face some missing libraries issues.

#### 1. The libreadline library is missing

 If you are getting the following error message while running APBS:

 ```
 apbs: error while loading shared libraries: libreadline.so.4: cannot open shared object file: No such file or directory 
 ```

 Your system is likely to have a `libreadline.so` in a more recent version. You can find it on a Linux distribution with the following command: `ldconfig -p | grep libreadline`

 On my machine, it tells me that `libreadline.so.8` is present in `/lib/x86_64-linux-gnu/`:
<pre><font color="#CC0000"><b>libreadline</b></font>.so.8 (libc6,x86-64) =&gt; /lib/x86_64-linux-gnu/<font color="#CC0000"><b>libreadline</b></font>.so.8
</pre>

So to fix the problem, i'll just have create a symlink to make APBS recognize `libreadline.so.8` as `libreadline.so.4`:

 ```
 sudo ln -s /lib/x86_64-linux-gnu/libreadline.so.8 /lib/x86_64-linux-gnu/libreadline.so.4
 ```

 #### 2. The libg2c library is missing

If you are getting the following error message while running APBS:

 ```
apbs: error while loading shared libraries: libg2c.so.0: cannot open shared object file: No such file or directory
 ```

 In this case, you'll have to install a gfortran compiler, for instance [g77](https://gcc.gnu.org/fortran/)

<br>
<br>

# How to cite SURFMAP
[Go to the top](#Table-of-contents)

If SURFMAP has been useful to your research, please cite us as well as the original MSMS paper:

> Hugo Schweke, Marie-Hélène Mucchielli, Nicolas Chevrollier, Simon Gosset, Anne Lopes (2021) SURFMAP: a software for mapping in two dimensions protein surface features. bioRxiv 2021.10.15.464543; doi: https://doi.org/10.1101/2021.10.15.464543

> Sanner MF, Olson AJ, Spehner JC. Reduced surface: an efficient way to compute molecular surfaces. Biopolymers. 1996 Mar;38(3):305-20. doi: 10.1002/(SICI)1097-0282(199603)38:3%3C305::AID-BIP4%3E3.0.CO;2-Y. PMID: 8906967. https://doi.org/10.1002/(sici)1097-0282(199603)38:3%3c305::aid-bip4%3e3.0.co;2-y

<br>

Moreover, if you use APBS in your research, please cite one or more of the following papers listed in the [Supporting APBS](https://apbs.readthedocs.io/en/latest/supporting.html) documentation page.
<br>
<br>

# References
[Go to the top](#Table-of-contents)

(1) Michel Sanner, Arthur J. Olson, Jean Claude Spehner (1996). Reduced Surface: an Efficient Way to Compute Molecular Surfaces. Biopolymers, Vol 38, (3), 305-320.

(2) Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112-128, 2018.
