---
title: "README installation SURFMAP"
author: "Hugo Schweke"
date: "February 8, 2021"
output: pdf_document
---

# Requirements

SURFMAP was developed under Linux on a Ubuntu 18.04 distribution. It should work under Mac OSX but has never been tested. It does not work under Windows.
It is written in python 2.7, R 3.6 and bash. It calls the MSMS (1) and APBS (2) softwares. The first is mandatory while the second is needed only to do electrostatics calculations.



# Installation
After downloading SURFMAP, you need to do the following:


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

2 - For the installation, please follow the instruction in the following page:

https://ssbio.readthedocs.io/en/latest/instructions/msms.html

3 - Once MSMS is installed, export the path to MSMS in your bashrc. For this do the following:

- open the `.bashrc` file in your favorite text editor and add the following line:

```bash
export MSMS=/path/to/MSMS/
```

Then type source `~/.bashrc` in the terminal. From now on if you type `echo $MSMS` in your terminal the output should be the path to MSMS.

4- There is a known bug in the file pdb_to_xyzr. You must do the following for it to function correctly:

- open the file with your favorite text editor.

- at line 31, replace the content `{nawk 'BEGIN {` with `{nawk -v pathmsms="$MSMS" 'BEGIN {`

- at line 34, replace the content `numfile = "./atmtypenumbers"` with `numfile = pathmsms"./atmtypenumbers"`

**III-a - Install APBS**

If you want to compute electrostatics potential, you will also need to install APBS. We recommend to install the pre-compiled binaries from the version 3.0.0 that can be found [here](https://github.com/Electrostatics/apbs/releases/tag/v3.0.0).

The install documentation from the pre-compiled binaries is accessible [there](https://apbs.readthedocs.io/en/latest/getting/index.html#installing-from-pre-compiled-binaries).


After downloading the pre-compiled binaries, you will need to edit your `~/.bashrc` file to make your system aware of the location of some APBS-3.0.0 required path. Let's say you have downloaded the (Linux) pre-compiled binaries of APBS here: `~/APBS-3.0.0.Linux`. Then open your favorite text editor and add the following lines:

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


## Docker command
```
docker run -it --rm -v $(pwd):/input/:ro -v $(pwd):/output/:rw surfmap_v0 -pdb /input/1vew_A.pdb -tomap electrostatics -d /output/output_elec
```

```
docker run -it --rm -v $(pwd):/home/surfmap/data/:rw surfmap_v0 -pdb /home/surfmap/data/1vew_A.pdb -tomap bfactor -d /home/surfmap/data/test
```

```
docker run -it --rm --entrypoint=/bin/bash surfmap_v0
```
