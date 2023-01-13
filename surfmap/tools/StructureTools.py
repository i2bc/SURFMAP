#!/usr/bin/env python3


# Anne Lopes
# module with several structure or PDB related functions...

import os, math

#########################################
#
#         deal PDB
#
###########################################

# parse PDB
def PDB_parser(infile, returnatomlist=False):
    """parses a pdb file and returns:
            - a pdb dict (with the following structure, chain->residue->atoms
            - a list (if returnatomlis = True) containing 3d coords of all atoms (useful for calcuating the CV)
            input: filename
            output : dPDB (PDB dict ), listXYZ (list containing 3d coords of all atoms)

            """
    # read pdb
    # ---------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var ini
    # ---------
    dPDB = {}
    atomlist = []
    dPDB["chains"] = []

    # loops over lines
    # -----------------
    for line in lines:
        if line[0:4] == "ATOM":

            # gets chain
            chain = line[21]

            # if not chain, creates the a new subdict for the chain and stores the chain in chainlist
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)  # add "chain" to chainlist
                dPDB[chain] = {}  # creates dict for the current chain
                # init list of residues for this chain
                dPDB[chain]["reslist"] = []

            # gets residue info
            curres = "%s" % (line[22:26]).strip()

            # if not res, creates the res key and adds the residue in reslist
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                # init atomlist for this res
                dPDB[chain][curres]["atomlist"] = []
                # gets resname
                dPDB[chain][curres]["resname"] = line[17:20].strip()

            # gets atom info
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
            dPDB[chain][curres][atomtype]["bfactor"] = line[60:67].strip()

            # modifs: creates a list and adds coords 3D of atome i in atomlist
            atomlist.append((float(line[30:38]), float(line[38:46]),
                             float(line[46:54])))  

    if returnatomlist == True:

        return dPDB, atomlist
    else:
        return dPDB


# parse PDB


def PDB_parser_dimdim(infile, reslimit, returnatomlist=False):
    """VERY SPECIFIC CASE FOR PARSING RESULTS FROM DOCKING WITH ATTRACT (chain A vs chain B)
        INVOLVING A FUSED DIMER (A) VS ANOTHER FUSED DIMER (B)
        parses a pdb file and returns:
            - a pdb dict (with the following structure, chain->residue->atoms
            - a list (if returnatomlis = True) containing 3d coords of all atoms (useful for calcuating the CV)
            input: filename
            output : dPDB (PDB dict ), listXYZ (list containing 3d coords of all atoms)

            """
    # read pdb
    # ---------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    # var ini
    # ---------
    dPDB = {}
    atomlist = []
    dPDB["chains"] = []

    # loops over lines
    # -----------------
    for line in lines:
        if line[0:4] == "ATOM":

            # gets chain
            chain = line[21]

            if chain == "A" :
                rec = True
            elif chain == "B" :
                rec = False
                
            # gets residue info
            curres = "%s" % (line[22:26]).strip()

            if int(curres) > reslimit and rec : # means we are in the next chain of the REC
                chain = "C"
            elif int(curres) > reslimit :
                chain = "D"
            else :
                chain = chain

            # if not chain, creates the a new subdict for the chain and stores the chain in chainlist
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)  # add "chain" to chainlist
                dPDB[chain] = {}  # creates dict for the current chain
                # init list of residues for this chain
                dPDB[chain]["reslist"] = []


            # if not res, creates the res key and adds the residue in reslist
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                # init atomlist for this res
                dPDB[chain][curres]["atomlist"] = []
                # gets resname
                dPDB[chain][curres]["resname"] = line[17:20].strip()

            # gets atom info
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
            dPDB[chain][curres][atomtype]["bfactor"] = line[60:67].strip()

            # modifs: creates a list and adds coords 3D of atome i in atomlist
            atomlist.append((float(line[30:38]), float(line[38:46]),
                             float(line[46:54])))  

    if returnatomlist == True:

        return dPDB, atomlist
    else:
        return dPDB


def ExtractSeqFromPDB(PDBfile) :
    """From a PDB file, returns a dico for which each key corresponds to 
       the chain and key(chain) points toward its corresponding sequence
       in one letter code. """

    f = open(PDBfile) 
    lines = f.readlines()
    f.close()

    # dict conatining the seq of each chain
    d_seq = {}

    # dict for 3 letter code to 1 letter code
    d_3to1 = {"ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E", "PHE" : "F",
              "GLY" : "G", "HIS" : "H", "ILE" : "I", "LYS" : "K", "LEU" : "L",
              "MET" : "M", "ASN" : "N", "PRO" : "P", "GLN" : "Q", "ARG" : "R",
              "SER" : "S", "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y"}
 
    # stores PDB in a dict
    dPDB = PDB_parser(PDBfile)

    # loops over chains
    for chaini in dPDB["chains"] :
        d_seq[chaini] = ""
        #loops over residues
        for resi in dPDB[chaini]["reslist"] :
            d_seq[chaini] += d_3to1[dPDB[chaini][resi]["resname"]]

    return d_seq

def initBfactor(dPDB, val = 0):
    """purpose: initiation of the bfactor key for each residue
       input: a dico with the dPDB format
    """

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"]:
            for atom in dPDB[chain][res]["atomlist"] :
                dPDB[chain][res][atom]["bfactor"] = val


#################################################
#           WRITING TOOLS
#################################################

def writePDB(dPDB, filout = "out.pdb", bfactor = False, chains = "all") :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each residue (one key per residue) in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    if chains == "all" :
        chainlist = dPDB["chains"]
    else:
        chainlist = [chains]

    for chain in chainlist:
        for res in dPDB[chain]["reslist"] :
            for atom in dPDB[chain][res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(
                        dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,
                        dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],
                        dPDB[chain][res][atom]["z"],dPDB[chain][res][atom]["bfactor"] ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(
                        dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],
                        chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],
                        dPDB[chain][res][atom]["z"] ))
                    
    fout.close()





def PDB2Fasta(PDBfile, outname = "no", splitfile = "no", outdir = "no") :
    """From a PDB file, returns a fasta file where each chain is translate into
       a one letter code sequence. 
       If splitfile == "yes", a fasta file will be generated for each chain, else
       a unique fasta file will be created containing all the sequences together.
       For example, A PDB file with 3 chains will lead to a 
       fasta file with three sequences corresponding to each chain."""


    d_seq = ExtractSeqFromPDB(PDBfile)
    nbchains = len(d_seq.keys())
    pdb = os.path.basename(PDBfile).split(".")[0]
    outlist = []

    if outdir == "no" :
        outdir = "outfasta"
        os.system("mkdir %s"%(outdir))

    if outname == "yes" and splitfile == "yes" :
        print ("These two options are not compatible, splitfile must be put off if outname is on.")
        return 0


    if outname == "no" :
        if splitfile == "no" :
            outlist.append[open("%s/%s.fasta"%(outdir, pdb), "w")] 
             
        else:
            for chain in d_seq :
                outlist.append(open("%s/%s%s.fasta"%(outdir, pdb,chain), "w"))
    else :
        outlist.append(open("%s"%outname, "w"))

    # loops over each chain and writes the >PDB id&chain and the corresponding seq
    for i in range(len(d_seq.keys())) :
        chain = list(d_seq.keys())[i]
        outlist[i].write(">%s:%s|PDBID|CHAIN|SEQUENCE\n"%(pdb[0:4], chain))
        for aa in d_seq[chain] :
            outlist[i].write("%s"%(aa))

        outlist[i].write("\n")

    # closes fasta files
    for out in outlist :
        out.close()



#####################################################
#         3D MANIPULATION TOOLS
#####################################################

def centerMassResidueList(dPDB, all = True, reslist = False, computeForMulti = True, chain = " "):
    """Calculates the center of mass of all the atoms contained in dPDB (all = True & reslist = False) or 
       for the atoms from a subset of residues given in the residue list (["12_A", "13_A", "27_A"])"""


    # check if PDB is multichains and if the user want to compute the CM considering all the chains
    #print ("nb chains ", len(dPDB["chains"]))
    #print ("computeForMulti ", computeForMulti)
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
            print ("please enter a reslist or fix all = True")
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



def spherical2CartCoord(ori, R, phi, theta):
    """purpose: from the coords of a center considered as the origin (ori) (tupple), the R distance between the 2 points
       (ori and atom for which we want to compute the cartesian coords)
       and a couple of angles phi, theta defining the position of the atom of interest,
       returns the coords of this atom (cartesian coords).
       input: ori (tupple corresponding to the 3D coords of the origin from which the spherical coords of
       the atom of interest have been defined)
              R distance between ori and the atom of interest
              phi, theta, angles defining the position of the atom of interest according to ori (spherical coordinates)
       output: cartcoords, a tupple containing the 3D cartesian coords xi, yi, zi of the atom of interest       
    """
    #print( "ori ", ori)
    #print ("R ", R)
    #print( "phi thet ", phi, theta)
    # computing x,y,z the 3D coords (relative to ori)
    z = math.cos(theta)*R
    x = math.sin(theta)*math.cos(phi)*R
    y = math.sin(phi)*math.sin(theta)*R

    # back to the original referential, storing in cartcoords
    cartcoords = (ori[0]+x, ori[1]+y, ori[2]+z) 
    
    return  cartcoords



def coord2spherical(ori, coordi):
    """purpose: from the 3D cartesian coordinates of an atom of interest (coordi) and
       those of a point considered as the origin (ori), computes its corresponding
       spherical coordinates according to the origin.
       input: ori (tupple corresponding to the cartesian coords of the origin)
              coordi (tupple corresponding to the cartesian coords of the input atom)
       output: spherical coords (R, phi, theta) of the input atom relative to ori       
    """

    #(x,y,z) = centering(ori, coordi)

    x = coordi[0] - ori[0]
    y = coordi[1] - ori[1]
    z = coordi[2] - ori[2]

    R = math.sqrt(x*x + y*y + z*z)
    theta = math.acos(z/R)

    # phi is given in [0:2pi]
    phi = 2*math.pi + math.atan2(y,x)
    phi = phi % (2*math.pi)

    return R, phi, theta



def computeDistance(d_res1, d_res2, mode="atom"):
    """Purpose: computes the distance between 2 residues (res1 and res2)
                Distance can be calculated either as the min distance between all the pairwise distances (atom-atom)
                or between the 2 centroids of the residues

       Inputs:  d_res1, d_res2 which are dico corresponding to res 1 and res 2 respectively
                mode: either "atom" or "centroid" (string)

       Output:  distance (float) """

    # mode atom: min distance atom-atom
    # --------------------------------------
    if mode == "atom":

        # *** computes all pairwise distances between atoms of res1, res2
        distances = []
        # loop over all atoms of res1
        for atom1 in d_res1["atomlist"]:
            coord1 = [d_res1[atom1]["x"], d_res1[atom1]["y"], d_res1[atom1]["z"]]
            # loop over all atoms of res2
            for atom2 in d_res2["atomlist"]:
                coord2 = [d_res2[atom2]["x"], d_res2[atom2]["y"], d_res2[atom2]["z"]]
                distances.append(distancePoints((coord1[0], coord1[1], coord1[2]), (coord2[0], coord2[1], coord2[2])))
        d_res1_res2 = min(distances)

    # mode centroid: distance between the centroids of the 2 given residues
    # ----------------------------------------------------------------------
    elif mode == "centroid":
        cent1 = getCentroid(d_res1)
        cent2 = getCentroid(d_res2)
        # print(cent1, cent2)

        d_res1_res2 = distancePoints(cent1, cent2)

    return d_res1_res2



def extractContactResidues(matdist, seuil) :
    """from a distance matrix (matdist), returns pairs of residues in contacts (seuil) in a list of lists """

    contacts = []
    for i in range(len(matdist[0])) :
        for j in range (i+1, len(matdist[0])) :
            if matdist[i][j] <= seuil :
                      contacts.append([i, j])

    return contacts


def distancePoints(L1, L2):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates
       output: distance"""

    # print(x1, x2)

    x = L1[0] - L2[0]
    y = L1[1] - L2[1]
    z = L1[2] - L2[2]
    return math.sqrt(x * x + y * y + z * z)


def getCentroid(d_res):
    """Purpose: Calculates the center of mass of a residue
       Input: a dico residue
       Output: coords x, y, z of the centroid (tuple format)
    """

    x = y = z = 0.0

    # loop over all atoms
    for atom in d_res["atomlist"]:
        x += d_res[atom]["x"]
        y += d_res[atom]["y"]
        z += d_res[atom]["z"]

    Xcen = float(x) / len(d_res["atomlist"])
    Ycen = float(y) / len(d_res["atomlist"])
    Zcen = float(z) / len(d_res["atomlist"])

    return (Xcen, Ycen, Zcen)


###########################
#### debug functions
###########################

def generateFastPDB(x, y, z, res = "GEN", atomname = "X", atomid = 1, resid = 1, chain = " ", bfactor = "", inser = " "):


    dPDB = {}
    dPDB["chains"] = [chain]
    dPDB[chain] = {}
    dPDB[chain]["reslist"] = [resid]
    dPDB[chain][resid] = {}
    dPDB[chain][resid]["atomlist"] = [atomname]
    dPDB[chain][resid][atomname] = {}
    dPDB[chain][resid][atomname]["id"] = atomid
    dPDB[chain][resid]["resname"] = res
    dPDB[chain][resid]["resnum"] = res
    dPDB[chain][resid]["inser"] = inser
    dPDB[chain][resid][atomname]["x"] = x
    dPDB[chain][resid][atomname]["y"] = y
    dPDB[chain][resid][atomname]["z"] = z
    if bfactor != "":
        dPDB[chain][resid][atomname]["bfactor"] = bfactor

    return dPDB

    
    
