# Overview

SURFMAP is a small software designed to compute different surface properties (see below) of a protein, and to map each of them on a 2-D plan through a sinusoidal projection. It enables the rapid visualization of a surface property across the whole surface of a protein. Such view can prove to be very useful when comparing the surface properties of homologous proteins, for example. 

SURFMAP was developed under Linux on a Ubuntu 18.04 distribution. It is written in python (version 3.7), R (version 3.6) and bash. It calls the MSMS (1) and APBS (2) softwares. The first is mandatory while the second is only required if the user wants to perform electrostatics calculations.

## Surface properties mapped by SURFMAP

SURFMAP provides the user the ability to compute 5 different protein surface properties: 

| Protein surface property | Description |
|----|----|
| Stickiness (1) | Measure of the propensity of an amino acid to be enriched (or depleted) at protein binding sites. The stickiness enables to detect regions that are theoretically more prone to interaction, i.e. "more sticky", from other regions. |
| Kyte-Doolittle hydrophobicity (2) | Measure of the degree of hydrophobicity/hydrophilicity of amino acids according to the Kyte-Doolittle hydrophobicity scale. |
| Wimley–White hydrophobicity (3) | Measure of the degree of hydrophobicity/hydrophilicity of amino acids according to the Wimley–White hydrophobicity scale. |
| Circular variance (4) | This measure characterizes the geometric properties of molecular structures by distinguishing between atoms accessible to the surface and buried atoms. In our case it is useful for the mapping of geometrical properties of the protein surface such as cavities and protuberant regions. |
| Electrostatic potential through APBS (5) | SURFMAP can call the APBS software that will compute the electrostatic potential, using the CHARMM forcefield (6), and map the resulting potential on a 2-D map. |
| B-factor values | The user have the possibility to map any value contained in the b-factor column of the pdb file provided in input. Two options are available: the mapping of either continuous values or discrete values. |



# Download

You can click [here](https://github.com/i2bc/SURFMAP/releases/latest/) to access the latest release.

After downloading, you'll need to uncompress the project. If you downloaded:
- the .zip file: type on a terminal (Unix and macOS): `unzip SURFMAP-x.x.x.zip`
- the .tar.gz file: type on a terminal (Unix and macOS): `tar xzvf SURFMAP-x.x.x.tar.gz`


Alternatively, you can also clone the whole project:
```
git clone https://github.com/i2bc/ORFmine.git
```


# Installation

There are two ways to install SURFMAP:
- the simplest and recommended way is to use a docker image of SURFMAP. It works on any system that can run docker (Unix, MacOS, or Windows 10 through WSL2).
- the alternative is to manually install SURFMAP on your system (only for Unix and MacOS). You will need to install MSMS (1) and APBS (2) (facultative) separately from SURFMAP.


## Use of SURFMAP through a docker image
You’ll first need to create an account on [docker hub](https://hub.docker.com/) and [install docker](https://docs.docker.com/get-docker/) on your machine. 

Once you’ve successfully registered to the hub and installed docker on your machine, you should be ready to use the docker image of SURFMAP through the script `run_surfmap.py`. First make sure this script is executable (the following command should also work on a Windows 10 machine through WSL2):
``` bash
chmod +x run_surfmap.py
```

Additionally you can create an alias of this python script to make it accessible from anywhere on your machine. To do so, add the following lines at the end of your `~/.bashrc` (or `~/.bash_profile` or `~/.profile`) file:

```bash
# Alias to run surfmap from its docker image"
alias surfmap='python3 PATH/TO/SURFMAP/run_surfmap.py'
```

Make sure to replace `PATH/TO/SURFMAP/` with the absolute path you’ve downloaded SURFMAP to.

From now on, you should be able to use SURFMAP by simply typing in a terminal `surfmap`. Note that the first time you’ll type it you may have to type docker login in the terminal and fill the fields with your docker ID and docker password. It will then take few minutes to download the image on your machine. 

Once the SURFMAP image has been successfully downloaded on your machine, you should see after typing the `surfmap` command in a terminal the following message:

<pre>
<font color="#4E9A06"><b>tutor@surfmap</b></font>:<font color="#3465A4"><b>~/i2bc/SURFMAP</b></font>$ surfmap
usage: run_surfmap.py [-h] -pdb PDB -tomap
                      {electrostatics,all,circular_variance,circular_variance_atom,wimley_white,stickiness,kyte_doolittle,binding_sites,bfactor}
                      [-coords COORDS] [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth]
                      [--png] [--keep]
run_surfmap.py: error: the following arguments are required: -pdb, -tomap
</pre>


## Manual installation of SURFMAP

**I - R packages and python library**

Install the following R packages inside an R console:

```R
install.packages("optparse")
install.packages("data.table")
install.packages("gplots")
```

Install the python library numpy, in Ubuntu you can do it with the following command:

```bash
sudo apt install python-numpy
```

**II - Install MSMS**

1 - To download the software, please go to the page:

http://mgltools.scripps.edu/downloads#msms

2 - For the installation, please follow the instructions in the following page:

https://ssbio.readthedocs.io/en/latest/instructions/msms.html

3 - Once MSMS is installed, export the path to MSMS in your bashrc. For this do the following:

- open the `.bashrc` file in your favorite text editor and add the following line:

```bash
export MSMS=/path/to/MSMS/
```

Then type `source ~/.bashrc` in the terminal to make the export effective.

From now on if you type the `msms` command in your terminal the output should be as shown below:

<pre><font color="#4E9A06"><b>tutor@surfmap</b></font>:<font color="#3465A4"><b>~/i2bc/SURFMAP</b></font>$ msms
MSMS 2.6.1 started on nchenche
Copyright M.F. Sanner (1994)
Compilation flags -O2 -DVERBOSE -DTIMING
MSMS: No input stream specified
</pre>


4- There is a known bug in the file pdb_to_xyzr. You must do the following for it to function correctly:

- open the file with your favorite text editor.

- at line 31, replace the content `{nawk 'BEGIN {` with `{nawk -v pathmsms="$MSMS" 'BEGIN {`

- at line 34, replace the content `numfile = "./atmtypenumbers"` with `numfile = pathmsms"./atmtypenumbers"`

**III-a - Install APBS**

If you want to compute electrostatics potential, you will also need to install APBS. We recommend to install the pre-compiled binaries from the version 3.0.0 that can be found [here](https://github.com/Electrostatics/apbs/releases/tag/v3.0.0).

The install documentation from the pre-compiled binaries is accessible [there](https://apbs.readthedocs.io/en/latest/getting/index.html#installing-from-pre-compiled-binaries).


After downloading the pre-compiled binaries, you will need to edit your `~/.bashrc` file to make your system aware of the location of some APBS-3.0.0 required paths. Let's say you have downloaded the (Linux) pre-compiled binaries of APBS here: `~/APBS-3.0.0.Linux`. Then open your favorite text editor and add the following lines:

```bash
export PATH="$HOME/APBS-3.0.0.Linux/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/APBS-3.0.0.Linux/lib"
export APBS="$HOME/APBS-3.0.0.Linux/share/apbs"
```

Finally type `source ~/.bashrc` in the terminal. From now on the `apbs` command should be accessible.

You now have to fix a bug in a file used by APBS. To do so, simply type the following command:
```bash
sed -i '152s/nsmall/nsmall\[i\]/' $APBS/tools/manip/psize.py
```

**III-b - Install pdb2pqr**

pdb2pqr is required to make easier the generation of the APBS input file. pdb2pqr is a python package that can be easily installed through the following command:
```bash
pip install pdb2pqr
```

A supplementary module is also required to use pdb2pqr:
```bash
pip install requests
```



Now SURFMAP should be ready for use. Open the MANUAL for more details about the software.


# References

(1) Michel Sanner, Arthur J. Olson, Jean Claude Spehner (1996). Reduced Surface: an Efficient Way to Compute Molecular Surfaces. Biopolymers, Vol 38, (3), 305-320.

(2) Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112-128, 2018.
