## Installation


### 1. Overview

ORFmine is a package that consists of two independent tools ORFtrack and ORFold. 
Both these tools have been developed in python3 (version >= 3.6).
The install.sh  script will install both ORFtrack and ORFold with their dependancies.
They can be used together or independently. 


### 2. Download and uncompress the latest release archive

##### Download the latest release
Here: [ ![](img/icons/download_16x16.png "Click to download the latest release")](https://github.com/i2bc/ORFmine/releases/latest/)
<br> or simply clone our github repository called [ORFmine](https://github.com/i2bc/ORFmine/)

##### Uncompress the archive
If you downloaded:

* the *.zip* file: ```unzip orfmine-x.x.x.zip```
* the *.tar.gz* file: ```tar xzvf orfmine-x.x.x.tar.gz```


### 3. Create an isolated environment
Although not strictly necessary, this step is highly recommended 
(it will allow you to work on different projects without having any conflicting library versions).
If you do not want to create a virtual environment, please go directly to the [install section](#general_install).
 
#### Install virtualenv
``` python
python3 -m pip install virtualenv
```

#### Create a virtual environment
```bash
virtualenv -p python3 orfmine_env
```

#### Activate the created environment
```bash
source orfmine_env/bin/activate
```

Once activated, any python library you will install using pip 
will be installed solely in this isolated environment.
You must activate this environment any time you need libraries installed 
in this environment. 

Once you are done working on your project, 
simply type `deactivate` to exit the environment.


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        To delete definitively your virtual environment, you can simply
        remove the directory with the following instruction:
        <code>rm -r my_env/</code>
    </p>
</div>

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
        We remind to the user that some external packages used in ORFmine 
	(such as Biopython) require python version >= 3.6. Before creating 
	your virtual environment make sure that your python version is up-to-date. 
    </p>
</div>

<a name="general_install"></a>

### 4. Install ORFMine 

#### Preparation before the Installation

If you just want to use **ORFtrack** in order to annotate all
the possible ORFs of a genome, you have no other dependencies 
to install, and you simply have to **Launch the Installation** 
presented [below](#launch_install). 

The installation of **ORFold** becomes a bit more demanding as
there are some external tools to be downloaded and/or installed 
before launching the installation.

Firstly, **ORFold** is based on the HCA method for the calculation of the
fold potential. As a result [pyHCA](https://github.com/T-B-F/pyHCA) 
[[1](https://www.biorxiv.org/content/10.1101/249995v1)]
is essential to be pre-installed in your machine before installing 
**ORFold**. You can [download](https://github.com/T-B-F/pyHCA)  for free and install **pyHCA** using 
the instructions of the developers.  
<br>
If you are not interested in the calculation of the disorder
and/or aggregation propensities with **ORFold** and you already
have installed pyHCA, you can simply launch the installation
presented [below](#launch_install).

However, in the case you want to use [IUPred](https://iupred2a.elte.hu) 
[2][3][4] and/or [Tango](http://tango.crg.es) [5][6][7] with **ORFold** you have to 
first contact their developers through the respective links and have access 
to their programs. These two softwares are not freely available for 
non-academic users.

Once you have access to the IUPred and Tango you have to place them in a directory
called ```softwares``` placed in the path: ```ORFmine/orfold_v1/orfold/```. To do so:


* First create the ```softwares``` directory if not already created:

```bash
mkdir ORFmine/orfold_v1/orfold/softwares
```

* Move the IUPred source code and data (provided by the developer):
	
		mv iupred2a.py ORFmine/orfold_v1/orfold/softwares
		mv data ORFmine/orfold_v1/orfold/softwares
	
* Move Tango source code:
	* For MacOS:
		
			mv tango2_3_1 ORFmine/orfold_v1/orfold/softwares

	* For linux:

			mv tango_x86_64_release ORFmine/orfold_v1/orfold/softwares

	* For windows:
		
			mv Tango.exe ORFmine/orfold_v1/orfold/softwares

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
    The calculation of the disorder or aggregation propensities are both optional and 
	complementary to the HCA score. As a result, IUPred and 
	Tango tools are not mandatory for the installation of ORFold. In addition,
	they are not necessarily coupled together. ORFold will properly be 
	installed without them or even with only one of them.    
    </p>
</div>
<a name="launch_install"></a>


#### Installation

If you use a virtual environment, be sure that your virtual environment is activated.
Then, in any case, follow the procedure described below:

 
```bash
cd ORFmine
chmod u+x install.sh
./install.sh
```

This script will first uninstall ORFmine if it was already installed and will
re-install it. In addition, it will install all the dependency packages needed for 
ORFtrack and ORFold.   

<br><br><br>
#### References

1. Bitard-Feildel, T. & Callebaut, I. HCAtk and pyHCA: A Toolkit and Python API for the Hydrophobic Cluster Analysis of Protein Sequences. bioRxiv 249995 (2018).
2. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
3. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
4. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
5. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
6. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004). 
7. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).

