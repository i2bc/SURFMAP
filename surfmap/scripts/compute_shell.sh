#!/bin/bash

# creator: Hugo Schweke
# email: hschweke@gmail.com
# date: 09/05/2021

# Script to compute a shell with MSMS from a pdb file.
# If elec option activated, call APBS to compute electrostatic potential at the specific positions around the surface calculated by MSMS.

set -e

# Absolute path to script
DIR="$(cd "$(dirname "$0")" && pwd)"

#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

# Set MSMS path
MSMS=$(dirname ${DIR})/utils/MSMS

#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`


function HELP {
    echo -e "\n**************************************************"
    echo -e "help documentation for ${SCRIPT}\n"
    echo -e "Script to compute apbs file from a simple pdb file.\n"
    echo "The following arguments are necessary:"
    echo "${NORM}${BOLD}-p --path input pdb file (file + path)"
    echo "-e --elec if 1 invoke APBS to calculate electrostatic potential, if 0 create simple shell."
    echo "-r --radius radius added to the standard radius of residues (the higher the radius, the smoother the surface, r=0 being the excluded surface area)."
    echo "-o --out main output directory."
    echo "-q --pqr Path to a PQR file"
    echo -e "-h --help${NORM}\n"
    exit 1
}

while getopts p:e:r:o:q:: option
do
    case "${option}" in
        p) pdb=${OPTARG};;
        e) elec=${OPTARG};;
        r) rad=${OPTARG};;
        o) outdir=${OPTARG};;
        q) pqr=${OPTARG};;
        h) HELP;;
        \?) # unrecognized option
        echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
        HELP
    esac
done


# Check number of arguments. If none are passed, print help and exit.
NUMARGS=$#

if [ $NUMARGS -eq 0 ]; then
    HELP
elif [ $NUMARGS -lt 6 ]; then
    echo -e "Missing argument(s) or syntax error -> showing help documentation.\n\n"
    HELP
fi


# pqr and apbs files are pdb file name with extension modified (respectively .pqr and .in files)
dshell=${outdir}/shells/
delec=${outdir}/tmp-elec/
pdbfile=$(basename $pdb)
xyzrfile=$dshell${pdbfile%.pdb}.xyzr # input needed for MSMS
vertfile=$dshell${pdbfile%.pdb}.vert # output of MSMS
csvfile=$dshell/${pdbfile%.pdb}.csv # input needed for multivalue
shellfile=$dshell/${pdbfile%.pdb}_shell.pdb # shell file with electrostatic values.

#============================= Shell creation =================================

# Create directory containing shells and other files
if [ ! -d "$dshell" ]
then
    mkdir -p $dshell
fi


# generate xyzr format input for MSMS
pdb2xyzr -pdb $pdb -rad $rad -out $xyzrfile

# Invoking MSMS to compute Connolly surface.
$MSMS/msms -if $xyzrfile -of $dshell/${pdbfile%.pdb} > log

# Modifying MSMS output file format (the .vert file that contains the coordinates of all the vertices of the surface model created).
# For this we just need to take the first three columns (x,y,z coordinates) and write them in a csv format file.
# This csv file is needed as an input for the multivalue executable (arguments = csv file and potdx file created by APBS).
tail -n +4 $vertfile | awk '{print $1","$2","$3}' > $csvfile


if [ "$elec" == 1 ]
then
    # create directory for electrostatics calculation
    if [ ! -d "$delec" ]
    then
        mkdir -p $delec
    fi

    # set io variables
    pqrfile=$delec/${pdbfile%.pdb}.pqr # pqr format is needed for APBS
    apbsfile=$delec/${pdbfile%.pdb}.in # input needed for APBS (in v2.0, changed pqr.in to .in)
    potfile=${pqrfile}.dx # output of APBS (OpenDX scalar format)
    multfile=$delec/${pdbfile%.pdb}.mult # output of multivalue (electrostatic potential calculated at the coordinates given in input).

    # generate or copy a given pqrfile in expected directory
    if [ -z "${pqr}" ]
    then
        pdb2pqr30 --ff CHARMM --whitespace $pdb $pqrfile >> log 2>&1 
    else 
        cp $pqr $delec/${pdbfile%.pdb}.pqr
    fi

    # Create apbs input file (.in)
    python3 $APBS/share/apbs/tools/manip/inputgen.py --method=auto $pqrfile > /dev/null 2>&1   # from APBS 3.4.1

    sed -i "s|${pdbfile%.pdb}.pqr|${pqrfile}|g" $apbsfile
    sed -i "s|pot dx pot|pot dx $pqrfile|g" $apbsfile

    # Invoke apbs -> create pot file
    $APBS/bin/apbs $apbsfile >> log 2>&1  # from APBS 3.4.1


    #=================== Electrostatic potential of the shell =====================

    # using "multivalue" executable from APBS to calculate electrostatic potential at the surface points determined by MSMS.
    # Multivalue is simply "crossing" the MSMS csv file with the potential gridfile created by APBS.
    $APBS/share/apbs/tools/bin/multivalue $csvfile $potfile $multfile >> log 2>&1  # from APBS 3.4.1
   
    # multivalue output converted into pdb format.
    multival_csv2pdb -i $multfile -o $shellfile

else
    multival_csv2pdb -i $csvfile -o $shellfile
fi

# Remove MSMS and APBS log files
if [ -f "log" ]; then
    rm log
fi
if [ -f "io.mc" ]; then
    rm io.mc
fi
