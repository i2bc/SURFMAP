#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, os.path, sys, string, re, math, numpy
from typing import Dict, List


PROPERTIES = {
    "stickiness": {"ALA":0.0062, "CYS":1.0372, "ASP":-0.7485, "GLU":-0.7893, "PHE":1.2727, "GLY":-0.1771, "HIS":0.1204, "ILE":1.1109, "LYS":-1.1806, "LEU":0.9138, "MET":1.0124, "ASN":-0.2693, "PRO":-0.1799 , "GLN":-0.4114, "ARG":-0.0876, "SER":0.1376, "THR":0.1031, "VAL":0.7599, "TRP":0.7925, "TYR":0.8806},
    "wimley_white": {"ALA":4.08, "CYS":4.49, "ASP":3.02, "GLU":2.23, "PHE":5.38, "GLY":4.24, "HIS":4.08, "ILE":4.52, "LYS":3.77, "LEU":4.81, "MET":4.48, "ASN":3.83, "PRO":3.80 , "GLN":3.67, "ARG":3.91, "SER":4.12, "THR":4.11, "VAL":4.18, "TRP":6.10, "TYR":5.19},
    "kyte_doolittle": {"ALA":1.8, "CYS":2.5, "ASP":-3.5, "GLU":-3.5, "PHE":2.8, "GLY":-0.4, "HIS":-3.2, "ILE":4.5, "LYS":-3.9, "LEU":3.8, "MET":1.9, "ASN":-3.5, "PRO":1.6 , "GLN":-3.5, "ARG":-4.5, "SER":-0.8, "THR":-0.7, "VAL":4.2, "TRP":-0.9, "TYR":-1.3}
}


def computeResidueCM(dPDB, chain, resnumber):
    ''' compute center of mass of a specific residue '''
    x = y = z = 0.0

    # looping over the current residue atoms
    for atom in dPDB[chain][resnumber]["atomlist"]:
        x +=dPDB[chain][resnumber][atom]["x"]
        y +=dPDB[chain][resnumber][atom]["y"]
        z +=dPDB[chain][resnumber][atom]["z"]
        
    Xcm = float(x)/len(dPDB[chain][resnumber]["atomlist"]) 
    Ycm = float(y)/len(dPDB[chain][resnumber]["atomlist"])
    Zcm = float(z)/len(dPDB[chain][resnumber]["atomlist"])

    return Xcm, Ycm, Zcm



def centerMassResidueList(dPDB, all = True, reslist = False, computeForMulti = True, chain = " "):
    """Calculates the center of mass of all the atoms contained in dPDB (all = True & reslist = False) or 
       for the atoms from a subset of residues given in the residue list (["12_A", "13_A", "27_A"])"""

    # check if PDB is multichains and if the user want to compute the CM considering all the chains
   # print "\ncenterMassResidueList"
   # print "nb chains ", len(dPDB["chains"])
   # print "computeForMulti ", computeForMulti, "\n"
    if (len(dPDB["chains"]) > 1) and (computeForMulti == True):
        compmulti = True
    elif (len(dPDB["chains"]) > 1) and (computeForMulti == False) :
        compmulti = False
        multi = True
    else:
        compmulti = False
        multi = False
        
    # case where we have SEVERAL chains and we want to compute the CM over ALL chains
    if compmulti == True : # we consider in that case ALL the residues of each chain and not a subset of res
 
       x = y = z = 0.0
       nbatoms = 0
        
       for chain in dPDB["chains"] :
           for res in dPDB[chain]["reslist"] :        

               # looping over the current residue atoms
               for atom in dPDB[chain][res]["atomlist"] :
                   x +=dPDB[chain][res][atom]["x"]
                   y +=dPDB[chain][res][atom]["y"]
                   z +=dPDB[chain][res][atom]["z"]
                   nbatoms +=1
            
       Xcm = float(x)/nbatoms
       Ycm = float(y)/nbatoms
       Zcm = float(z)/nbatoms
       return Xcm, Ycm, Zcm

    # case where we have ONE chain or SEVERAL chains but we want to compute CM over ONE chain given by arg 'chain'
    else :
        if multi == False : # means we have one chain
            chain = dPDB["chains"][0]
        # else it means we have several chains but we want to compute CM over one specific chain
            
        if all == True :
            reslist = dPDB[chain]["reslist"]
        elif (all != True) and (reslist == False):
            print("please enter a reslist or fix all = True")
            sys.exit()
            
        x = y = z = 0.0
        nbatoms = 0
        for res in reslist :        
            # looping over the current residue atoms
            for atom in dPDB[chain][res]["atomlist"] :
                x +=dPDB[chain][res][atom]["x"]
                y +=dPDB[chain][res][atom]["y"]
                z +=dPDB[chain][res][atom]["z"]
                nbatoms +=1

        Xcm = float(x)/nbatoms
        Ycm = float(y)/nbatoms
        Zcm = float(z)/nbatoms
        return Xcm, Ycm, Zcm



def get_coords_idres(dPDB: Dict):
    # create numpy array of coordinates (for speed purpose)
    coordlist = list()
    idres = list()
    for chaini in dPDB["chains"] :
        for resi in dPDB[chaini]["reslist"]:
            for atomi in dPDB[chaini][resi]["atomlist"] :
                xres = dPDB[chaini][resi][atomi]["x"]
                yres = dPDB[chaini][resi][atomi]["y"]
                zres = dPDB[chaini][resi][atomi]["z"]
                coordlist.append((xres,yres,zres))
                idres.append((chaini, resi, atomi))
    coordlist = numpy.asarray(coordlist)

    return coordlist, idres


def getAtomCMRDist(coordlist, idres, CMR) :
    """computes the distance between each atom stored in dPDB and the Center of Mass of given in input (CMR -> tupple)
       Returns a dico (Key = "Chaini_Resi_Atomi", Value = distance with CMR)
    """
    mindist = closest_atom_2(CMR, coordlist)
    resmini = "%s_%s_%s"%(idres[mindist])
    
    return mindist, resmini


def closest_atom(part, atoms):
    dist_2 = numpy.sum((atoms - part)**2, axis=1)
    return numpy.argmin(dist_2)


def closest_atom_2(part, atoms):
    deltas = atoms - part
    dist_2 = numpy.einsum('ij,ij->i', deltas, deltas)
    return numpy.argmin(dist_2)


################################################################
#               ENERGY, POTENTIAL AND PROPENSITY
################################################################


def CoulombicEner(qi, qj, rij, epsilon="default") :    
    if epsilon == "default" :
        epsilon = 15*rij

    return float(qi*qj)/(epsilon*rij)


def returnPropensity(aa, scale="stickiness"):    
    return PROPERTIES[scale][aa]

    
################################################################
#                DEALING WITH PDBs
################################################################


def parsePDBMultiChains(infile, charge = 1, chargeFromInfile = False, bfactor = False, CG = False) :
    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["reslist"] = []
    dPDB["chains"] = []
    mbf = 0
    
    # parcoure le PDB   
    for line in lines :
        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and ( line[17:20].strip() == "MET") or  ( line[17:20].strip() == "MSE") ) :
            chain = line[21]
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s"%(line[22:27]).strip()
            resnum = "%s"%(line[22:26]).strip()
            if not curres in dPDB[chain]["reslist"] : # first time we encounter it
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = line[17:20].strip()
                dPDB[chain][curres]["atomlist"] = []
                dPDB[chain][curres]["atomlistTowrite"] = []
                alternateoccupancy = None #"%s"%(line[16:17])
                occupancy = "%s"%(line[16:17]) 
                if occupancy != " " :
                    alternateoccupancy = occupancy
            else: # this is not a new residue
                occupancy = "%s"%(line[16:17])
                if occupancy != " " and alternateoccupancy == None : # means we are in the first alternate location of that residue
                    alternateoccupancy = occupancy
            if CG : # means we are parsing a CG model so we have to treat the CSE atomtypes which can be redondant in terms of name the same res
                atomtype = "%s_%s"%(line[6:11].strip(), line[12:16].strip())
            else:
                atomtype = line[12:16].strip()
            #if not atomtype in dPDB[chain][curres]["atomlist"] :
            if occupancy == alternateoccupancy  or occupancy == " " : # means this atom corresponds to the first rotamer found in the PDB for this residue
                if CG :
                    dPDB[chain][curres]["atomlistTowrite"].append(atomtype.split("_")[1]) # necessary for the writing later
                dPDB[chain][curres]["atomlist"].append(atomtype)
                dPDB[chain][curres][atomtype] = {}
                dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
                dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
                dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
                dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
                if bfactor:
                    try:
                        dPDB[chain][curres][atomtype]["bfactor"] = float(line[60:67].strip())
                    except ValueError:
                        dPDB[chain][curres][atomtype]["bfactor"] = 0
                        mbf += 1
                    #print "bfactor: ", float(line[60:67].strip())

                if chargeFromInfile == True :
                    dPDB[chain][curres][atomtype]["charge"] = float(line[60:67])
                else :
                    dPDB[chain][curres][atomtype]["charge"] = charge

            dPDB[chain][curres]["resnum"] = resnum
            dPDB[chain][curres]["inser"] = "%s"%(line[26:27])
            #print "cures ", curres
            #print dPDB[chain][curres]

    
    if mbf == 1:
        print("WARNING: There is " + str(mbf) + " missing bfactor in the input pdb:\nmissing bfactors are arbitrary set to zero\nThis is probably not what you want, you should check that the bfactor column is correctly set.\n")
    elif mbf > 1:
        print("WARNING: There are " + str(mbf) + " missing bfactors in the input pdb:\nmissing bfactors are arbitrary set to zero\nThis is probably not what you want, you should check that the bfactor column is correctly set.\n")
    
    return dPDB



def parsePDBParticule(filin, charge = 1, infile = False) :

    # lecture du fichier PDB 
    f = open(filin, "r")
    lines = f.readlines()
    f.close()

    # var init
    dPDB = {}
    dPDB["partlist"] = []
    
    # parcoure le PDB   
    for line in lines :
        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and ( (line[17:20].strip() == "MET") or  (line[17:20].strip() == "MSE") )) :
            numpart = line[7:11].strip()
            dPDB["partlist"].append(numpart)
            dPDB[numpart] = {}
            dPDB[numpart]["x"] = float(line[30:38])
            dPDB[numpart]["y"] = float(line[38:46])
            dPDB[numpart]["z"] = float(line[46:54])
            if infile == True :
                #print "line:" + line[61:67]
                dPDB[numpart]["charge"] = float(line[60:67])
            else :
                dPDB[numpart]["charge"] = charge

    return dPDB

    

def writePDB(dPDB, filout="out.pdb", bfactor=False, charge=False, CG=False, chains_to_rm: List=[]):
    """according to the coordinates in dPDB, writes the corresponding PDB file."""

    fout = open(filout, "w")
    #print dPDB["reslist"][1], dPDB[dPDB["reslist"][1]]["C"]["id"]

    for chain in dPDB["chains"]:
        if chain in chains_to_rm:
            continue
        
        for res in dPDB[chain]["reslist"] :
            
            if CG :
                listatomtowrite = dPDB[chain][res]["atomlistTowrite"]
            else: 
                listatomtowrite = dPDB[chain][res]["atomlist"] 

            listatomtoread = dPDB[chain][res]["atomlist"] 
            cmt = 0

            for atom in listatomtoread :
                
                if (bfactor == True) and (charge == True) :
                    print("You cannot store the charge and bfactor information. You have to choose!")
                    sys.exit()
                elif (bfactor == True) :
                    #print "BFACT", dPDB[chain][res][atom]["bfactor"]
                    fout.write("ATOM  %5s  %-4s%3s %s%4s%s   %8.3f%8.3f%8.3f  1.00%6.2f X X\n"%(dPDB[chain][res][atom]["id"], listatomtowrite[cmt], dPDB[chain][res]["resname"],chain, res,dPDB[chain][res]["inser"],dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"],dPDB[chain][res][atom]["bfactor"] ))

                elif (charge == True) :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s%s   %8.3f%8.3f%8.3f  1.00%6.2f X X\n"%(dPDB[chain][res][atom]["id"], listatomtowrite[cmt], dPDB[chain][res]["resname"],chain, res,dPDB[chain][res]["inser"],dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"],dPDB[chain][res][atom]["charge"] ))

                else:
                   fout.write("ATOM  %5s  %-4s%3s %s%4s%s   %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[chain][res][atom]["id"], listatomtowrite[cmt], dPDB[chain][res]["resname"],chain, res,dPDB[chain][res]["inser"],dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"] ))

                cmt+=1    
    fout.close()

         

############################################
#         3D MANIPULATIONS
############################################


def coord2spherical(ori, coordi):
    """purpose: from the 3D cartesian coordinates of an atom of interest (coordi) and
       those of a point considered as the origin (ori), computes its corresponding
       spherical coordinates according to the origin.
       input: ori (tupple corresponding to the cartesian coords of the origin)
              coordi (tupple corresponding to the cartesian coords of the input atom)
       output: spherical coords (R, phi, theta) of the input atom relative to ori       
    """
    x = coordi[0] - ori[0]
    y = coordi[1] - ori[1]
    z = coordi[2] - ori[2]

    R = math.sqrt(x*x + y*y + z*z)
    theta = math.acos(z/R)

    # phi is given in [0:2pi]
    phi = 2*math.pi + math.atan2(y,x)
    phi = phi % (2*math.pi)

    return R, phi, theta



###################### PART CV ##########################


#Compute distance between two atoms
def distance(atom1, atom2):
    x=atom1[0]-atom2[0]
    y=atom1[1]-atom2[1]
    z=atom1[2]-atom2[2]
    return math.sqrt(x**2+y**2+z**2)


#Compute list of atoms in a radius rc of the center of a reference atom
def Env_i(atomRef, rc, atoms):
    neighbor_atoms=[]
    for i in range(len(atoms)):
        if distance(atomRef, atoms[i])<rc:
            neighbor_atoms.append(atoms[i])
    return neighbor_atoms


#Compute vector defined by 2 atoms
def Rij(ati, atj):
    return (atj[0]-ati[0], atj[1]-ati[1], atj[2]-ati[2])


#Compute norm
def Norme(x,y,z):
    return math.sqrt(x**2+ y**2+ z**2)


#Compute CV of an atom
def CVi(ati, rc, atoms):
    neighbor_ats=Env_i(ati, rc, atoms)
    SVX=0.0
    SVY=0.0
    SVZ=0.0
    for j in range(len(neighbor_ats)):
        vx,vy,vz=Rij(ati, neighbor_ats[j])
        norm=Norme(vx,vy,vz)
        if norm!=0:
            SVX+=vx/norm
            SVY+=vy/norm
            SVZ+=vz/norm

    CV=1-Norme(SVX,SVY,SVZ)/len(neighbor_ats)
    if Norme(SVX,SVY,SVZ)/len(neighbor_ats) > 1:
        print("CV < 0 !")
        print("Norme: ", Norme(SVX,SVY,SVZ))
        print("nb close atoms: ", len(neighbor_ats),"\n")
    return CV


# Compute CV for all atoms of a list
def CV_AllRes(atoms, rc=8):
    allCVs=[]
    for i in range(len(atoms)):
        allCVs.append([atoms[i], CVi(atoms[i], rc, atoms)])
    return allCVs


# Compute circular variance of a pdb file
def compute_CV(pdb, perres=True, outfilename: str=None) :
    pdbfile = open(pdb, "r")
    pdblines = pdbfile.readlines()
    pdbfile.close()
    
    atomlist = []
    cpt = 0
    
    for line in pdblines:
        if line.startswith("ATOM") == True:
            cpt = cpt+1
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            res = line[17:20]
            numres = int(line[22:26])
            chain = line[21]
            atomlist.append((x, y, z, numres, res, chain))

    # Compute circular variance for all atom in file
    allCV=CV_AllRes(atomlist)
    
    dico_CV = dict()
    for elem in allCV:
        res = elem[0][4] + "_" + elem[0][5] + "_" + str(elem[0][3])
        CV = elem[1]

        if not res in dico_CV.keys():
            dico_CV[res] = list()
            dico_CV[res].append(list())
            dico_CV[res][0].append(CV)
        else:
            dico_CV[res][0].append(CV)

    for key in dico_CV:
        dico_CV[key].append(sum(dico_CV[key][0])/len(dico_CV[key][0]))


    if not outfilename:
        outfilename = os.path.splitext(pdb)[0] + "_CV.pdb"

    with open(outfilename, "w+") as outfile:
        cpt = 0
        for line in pdblines:
            if line.startswith("ATOM") == True:
                resid = line[17:20].strip(" ")+"_"+ line[21]+ "_" +line[22:26].strip(" ")
                toadd = 61-len(line)
                if toadd > 0:
                    white_spaces = " " * toadd
                else:
                    white_spaces = ""
                
                if perres == True:
                    CV = round(dico_CV[resid][1],2)
                else:
                    CV = round(allCV[cpt][1],2)
                    CV = '{0:.2f}'.format(CV)
                CV_to_write = '{:>6}'.format(CV)
                line_to_write = line[0:60].strip("\n") + white_spaces + CV_to_write + "\n"
                outfile.write(line_to_write)
                cpt = cpt+1
