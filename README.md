# Quick overview

SURFMAP is a free standalone and easy-to-use software that enables the fast and automated 2-D projection of either predefined features of protein surface (electrostatic potential, Kyte-Doolittle hydrophobicity, stickiness and surface relief) or any descriptor encoded in the temperature factor column of a PDB file. The 2-D maps computed by SURFMAP can thus be used to analyze and/or compare protein surface properties.


### Requirements

SURFMAP is written in python (version 3.7), R (version 3.6) and bash. It relies on the MSMS  software(1) and may optionally requires APBS (2) if the user wants to perform electrostatics calculations.

All these requirements (including the APBS software) are fullfilled in a [pre-built docker image of SURFMAP](#Run-SURFMAP-through-its-pre-built-docker-image) that we recommend the user to use.


# Download

You can click [here](https://github.com/i2bc/SURFMAP/releases/latest/) to access the latest release.

After downloading, you'll need to unarchive the project. If you downloaded:
- the .zip file: type on a terminal (Unix and macOS): `unzip SURFMAP-x.x.x.zip`
- the .tar.gz file: type on a terminal (Unix and macOS): `tar xzvf SURFMAP-x.x.x.tar.gz`


Alternatively, you can also clone the whole project:
```
git clone https://github.com/i2bc/ORFmine.git
```

After having unarchived (or cloned) the project, you should have the following files/directories:

<pre><font color="#3465A4"><b>.</b></font>
├── <font color="#3465A4"><b>example</b></font>
├── README.md
├── requirements.txt
├── <font color="#4E9A06"><b>run_surfmap.py</b></font>
├── <font color="#3465A4"><b>scripts</b></font>
├── <font color="#4E9A06"><b>SURFMAP_launcher.py</b></font>
├── SURFMAP_manual.pdf
├── SURFMAP_manual.Rmd
└── <font color="#3465A4"><b>tools</b></font>
</pre>


# Run SURFMAP through its pre-built docker image

The easiest and recommended way to run SURFMAP is to take profit of its pre-built docker image. This image will work on any system that can run docker (Unix, MacOS, or Windows 10 through WSL2). It includes all the dependencies and external softwares (MSMS & APBS) to run SURFMAP so that the user don't need to install anything on its machine (except docker).  

Yet if you want/need to install SURFMAP on your machine, please refer to the [fully manual installation](#Installation) section.  

### 1. Get docker on your machine

You’ll first need to create an account on [docker hub](https://hub.docker.com/) and [install docker](https://docs.docker.com/get-docker/) on your machine.


### 2. Make the proxy script `run_surfmap.py` executable from anywhere

Once you’ve successfully registered to the hub and installed docker on your machine, you should be ready to use the docker image of SURFMAP through the script `run_surfmap.py`.

First make sure this script is executable (the following command should also work on a Windows 10 machine through WSL2):
```bash
chmod +x run_surfmap.py
```

<a id="surfmap-alias"></a>
Additionally you can create an alias of this python script to make it accessible from anywhere on your machine. To do so, add the following lines at the end of your `~/.bashrc` (or `~/.bash_profile` or `~/.profile`) file:

```bash
# Alias to run surfmap from its docker image"
alias surfmap='python3 PATH/TO/SURFMAP/run_surfmap.py'
```

Make sure to replace `PATH/TO/SURFMAP/` with the absolute path you’ve downloaded SURFMAP to. Then type `source ~/.bashrc` (or `~/.bash_profile` or `~/.profile`) in the terminal to make the alias effective.

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

<br>


# Usage of SURFMAP

In the following section, we will assume that the SURFMAP program is called through the `surfmap` alias command (see [here](#surfmap-alias)). 

To illustrate how to use SURFMAP, we will use files in the `example/` directory that can be found in the downloaded SURFMAP project:

<pre><font color="#3465A4"><b>example/</b></font>
├── 1g3n_A_CV.pdb
├── 1g3n_A.pdb
├── 1gv3_A_binding_sites.pdb
├── <font color="#3465A4"><b>example_1g3n_A_stickiness</b></font>
├── <font color="#3465A4"><b>example_1g3n_A_stickiness_allfiles</b></font>
├── <font color="#3465A4"><b>example_1g3n_A_stickiness_mapping_residues</b></font>
├── <font color="#3465A4"><b>example_1gv3_A_stickiness_binding_sites</b></font>
├── <font color="#3465A4"><b>example_mapping_stickiness</b></font>
├── README
└── residues_to_map.txt
</pre>

### Mandatory arguments

SURFMAP requires 2 mandatory arguments as inputs:
<table style="margin-left: auto;  margin-right: auto;">
    <tr>
        <td><code>-pdb</code></td>
        <td>must be followed by a protein structure in a PDB format</td>
    </tr>
    <tr>
        <td><code>-tomap</code></td>
        <td>must be followed by one of the following options:
            <ul style="list-style-type:none">
                <li><code>wimley_white</code></li>
                <li><code>kyte_doolittle</code></li>
                <li><code>stickiness</code></li>
                <li><code>circular_variance</code></li>
                <li><code>circular_variance_atom</code></li>
                <li><code>electrostatics</code></li>
                <li><code>bfactor</code></li>
                <li><code>binding_sites</code></li>
                <li><code>all</code></li>
            </ul>
        </td>
    </tr>    
</table>

| |  |
| --- | --- |
| `-pdb` | must be followed by a protein structure in a PDB format |
| `-tomap` | must be followed by one of the following options:<br>&nbsp;&nbsp;- wimley_white<br>&nbsp;&nbsp;- kyte_doolittle<br>&nbsp;&nbsp;- stickiness<br>&nbsp;&nbsp;- circular_variance<br>&nbsp;&nbsp;- circular_variance_atom<br>&nbsp;&nbsp;- electrostatics<br>&nbsp;&nbsp;- bfactor<br>&nbsp;&nbsp;- binding_sites<br>&nbsp;&nbsp;- all |



# Fully manual installation of SURFMAP

### Requirements

- a Unix-based system machine (any Linux distribution, MacOS, or Windows 10 through WSL2)
- python >= 3.7
- R >= 3.6
- external softwares:
  - MSMS
  - APBS (only if you want to perform electrostatics calculations)

## I. Install the required R packages and python library

Please open an R console and install the following R packages:

```R
install.packages("optparse")
install.packages("data.table")
install.packages("gplots")
```

Install the python dependencies with the following command:

```bash
pip install -r requirements.txt
```

## II. Install MSMS

#### 1. Get the MSMS software

You can download MSMS [here](http://mgltools.scripps.edu/downloads#msms)


#### 2. Indicate your system where MSMS is located

After downloading MSMS and having unarchive it, open the `.bashrc` (or `~/.bash_profile` or `~/.profile`) file in your favorite text editor and add the following line:

```bash
export MSMS=/path/to/MSMS/
```

Then type `source ~/.bashrc` (or `~/.bash_profile` or `~/.profile`) in the terminal to make the export effective.

From now on if you type the `msms` command in your terminal the output should be as shown below:

<pre><font color="#4E9A06"><b>tutor@surfmap</b></font>:<font color="#3465A4"><b>~/i2bc/SURFMAP</b></font>$ msms
MSMS 2.6.1 started on nchenche
Copyright M.F. Sanner (1994)
Compilation flags -O2 -DVERBOSE -DTIMING
MSMS: No input stream specified
</pre>

Note: If anythong goes wrong, you can refer to [MSMS instruction page](https://ssbio.readthedocs.io/en/latest/instructions/msms.html).


#### 3. Correct a known bug in the pdb_to_xyzr file 

Please open `pdb_to_xyzr` the file with your favorite text editor and then:

- at line 31, replace the content `{nawk 'BEGIN {` with `{nawk -v pathmsms="$MSMS" 'BEGIN {`

- at line 34, replace the content `numfile = "./atmtypenumbers"` with `numfile = pathmsms"./atmtypenumbers"`

## III. Install APBS

If you want to compute electrostatics potential, you will also need to install APBS. We recommend to install the pre-compiled binaries from the version 3.0.0 that can be found [here](https://github.com/Electrostatics/apbs/releases/tag/v3.0.0).

The install documentation from the pre-compiled binaries is accessible [there](https://apbs.readthedocs.io/en/latest/getting/index.html#installing-from-pre-compiled-binaries).


After downloading the pre-compiled binaries, you will need to edit your `~/.bashrc` (or `~/.bash_profile` or `~/.profile`) file to make your system aware of the location of some APBS-3.0.0 required paths. Let's say you have downloaded the (Linux) pre-compiled binaries of APBS here: `~/APBS-3.0.0.Linux`. Then open your favorite text editor and add the following lines:

```bash
export PATH="$HOME/APBS-3.0.0.Linux/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/APBS-3.0.0.Linux/lib"
export APBS="$HOME/APBS-3.0.0.Linux/share/apbs"
```

Finally type `source ~/.bashrc` (or `~/.bash_profile` or `~/.profile`) in the terminal. From now on the `apbs` command should be accessible.

You now have to fix a bug in a file used by APBS. To do so, simply type the following command:
```bash
sed -i '152s/nsmall/nsmall\[i\]/' $APBS/tools/manip/psize.py
```

Now SURFMAP should be ready for use. After typing in a terminal `apbs` you shoul see:

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




# References

(1) Michel Sanner, Arthur J. Olson, Jean Claude Spehner (1996). Reduced Surface: an Efficient Way to Compute Molecular Surfaces. Biopolymers, Vol 38, (3), 305-320.

(2) Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Science, 27, 112-128, 2018.
