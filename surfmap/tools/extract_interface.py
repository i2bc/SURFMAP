#!/usr/bin/env python

# extract_interface.py
# purpose: extracts binding sites from a dimer
# Anne Lopes 
# must install freesasa with pip before
import os, sys
import StructureTools
#import CV, CV2
import freesasa
# Def functions
#==============
def usage():
    print ("""

     obligatory:
     ===========

     -pdb     -> pdb of the cplx


     optional:
     =========

     -mono   -> dir containing the pdb of the monomers (as exctracted from the cplx, ie, bound state, not unbound!)
                If calling this option, monomers must follow the following pattern
                 [rootname_cpl_pdb]_r.pdb and [rootname_cpl_pdb]_l.pdb for the rec and lig respectively
                (default: no monomers, they will be generated from the complex automatically)

    -bfact   -> writes the binding site residues in the bfactor of each mono (default: False)

     
     
  """)


def getSurfRes(dASAcp, dASAlig, dASArec):
    """extract surface residues through the comparison of their ASA in the monomer VS complex
    takes as inputs: dico ASA for cplx, lig and rec respectively

    """
        
    #init var
    ligsurf = []
    liginter = []
    recsurf = []
    recinter = []

    chainlig = list(dASAlig.keys())[0]
    chainrec = list(dASArec.keys())[0]

    # loop over cplx chains
    for chncp in dASAcp :

        # loop over residues                 
        for rescp in dASAcp[chncp] :
                
            if chncp == chainlig : #means the current chain is the lig
                        
                dASA = dASAlig[chncp][str(rescp)].relativeTotal - dASAcp[chncp][str(rescp)].relativeTotal
                    
                if dASA > 0 : # means the ASA in the lig is greater than in the cplx form, so the res is buried in the interface                        
                    liginter.append(rescp)
                    ligsurf.append(rescp)
                else: # does not belong to the inter
                    if dASAlig[chncp][str(rescp)].relativeTotal > 0.25 : # means the res is located on the lig surface
                        ligsurf.append(rescp)
            else : # dealing the rec
                dASA = dASArec[chncp][str(rescp)].relativeTotal - dASAcp[chncp][str(rescp)].relativeTotal
                    
                if dASA > 0 : # means the ASA in the lig is greater than in the cplx form, so the res is buried in the interface                        
                    recinter.append(rescp)
                    recsurf.append(rescp)
                else: # not in interface
                    if dASArec[chncp][str(rescp)].relativeTotal > 0.25 : # means the res is located on the lig surface
                        recsurf.append(rescp)              
                    
            
    return recsurf, ligsurf, recinter, liginter

def write_interRes_inBfactor(dPDBrec, chainrec, recinter, dPDBlig, chainlig, liginter):
    
    # stores status in bfactro for REC
    for resi in dPDBrec[chainrec]["reslist"] :
        if resi in recinter :
            bf = 1
        else:
            bf = 0
        for atomi in dPDBrec[chainrec][resi]["atomlist"] :
            dPDBrec[chainrec][resi][atomi]["bfactor"] = bf

    # stores status in bfactro for REC
    for resi in dPDBlig[chainlig]["reslist"] :
        if resi in liginter :
            bf = 1
        else:
            bf = 0
        for atomi in dPDBlig[chainlig][resi]["atomlist"] :
            dPDBlig[chainlig][resi][atomi]["bfactor"] = bf

    # writes corresponding pdb
    StructureTools.writePDB(dPDBrec, filout="%s/%s_r_inter_bf.pdb"%(monodir, rootname), bfactor=True)
    StructureTools.writePDB(dPDBlig, filout="%s/%s_l_inter_bf.pdb"%(monodir, rootname), bfactor=True)
    


# Get Arguments
#===============


try:
    cplxfile = sys.argv[sys.argv.index("-pdb")+1]
    #print ("cplx to treat:", cplxfile)
    #print("caution, deals with dimers only! Returns the binding sites involved in the input dimer")
    rootname = os.path.basename(cplxfile).split(".")[0]
except:    
    usage()
    print ("ERROR: please, enter the pdbfile of the cplx to treat")
    sys.exit()

try:
    monodir = sys.argv[sys.argv.index("-mono")+1]
    #print ("monomers to treat are stored in ", monodir)
    #print("caution, monomers must follows the name pattern [cplx_name]_l.pdb, [cplx_name]_r.pdb")
    generate = False
except:    
    generate = True

try:
    sys.argv.index("-bfact")
    #print ("binding site residues will be stored in the bfactor ")
    bfact = True
except:    
    bfact = False


try:
    sys.argv.index("-writeLenInterOnly")
    writeLenInterOnly = True
except:    
    writeLenInterOnly = False


# prepare data
dico_inter_surf = {}
dCPLX = StructureTools.PDB_parser(cplxfile)
if len(dCPLX["chains"]) != 2 :
    #print("Input error, deals with dimers only, please provide a dimer!")
    sys.exit()

# if generate == True, generates the pdb of each monomer, else will take those provided in -mono 
if generate:
    os.makedirs("monomers", exist_ok=True)
    monodir = "monomers"
    dPDBrec = {}
    chainrec = dCPLX["chains"][0]
    dPDBrec["chains"]= [chainrec]
    dPDBrec[chainrec] = dCPLX[chainrec]
    dPDBlig = {}
    chainlig = dCPLX["chains"][1]
    dPDBlig["chains"]= [chainlig]
    dPDBlig[chainlig] = dCPLX[chainlig]

    #if not writeLenInterOnly:
    # write the pdb for each monomers from the corresponding dict
    StructureTools.writePDB(dPDBrec, filout="%s/%s_r.pdb"%(monodir, rootname), bfactor=False)
    StructureTools.writePDB(dPDBlig, filout="%s/%s_l.pdb"%(monodir, rootname), bfactor=False)

# path and name of lig and rec (monomers)
monol = "%s/%s_l.pdb"%(monodir, rootname)       
monor = "%s/%s_r.pdb"%(monodir, rootname)       

#print("dealing with monomers %s and %s "%(monol, monor))
    
# init dico for this complex
dico_inter_surf["lig"] = {}
dico_inter_surf["rec"] = {}

# compute SASA for the cplx
strcplx = freesasa.Structure("%s"%(cplxfile))
outASAcplx = freesasa.calc(strcplx)
dASAcp = outASAcplx.residueAreas()

# compute SASA for the lig
strlig = freesasa.Structure("%s"%(monol))
outASAlig = freesasa.calc(strlig)
dASAlig = outASAlig.residueAreas()

# compute SASA for the rec
strrec = freesasa.Structure("%s"%(monor))
outASArec = freesasa.calc(strrec)
dASArec = outASArec.residueAreas()

if not generate : # must extract the names of each chain
    chainlig = list(dASAlig.keys())[0]
    chainrec = list(dASArec.keys())[0]

recsurf, ligsurf, recinter, liginter = getSurfRes(dASAcp, dASAlig, dASArec)


if writeLenInterOnly:
    out = open("leninter.txt", "w")
    out.write("%s %s\n"%(len(recinter), len(liginter)))
    out.close()
    sys.exit()
    

# writes the residues of the binding sites in a text file + bfactor if bfact == True
out_r = open("%s/%s_%s_bs.txt"%(os.path.dirname(monor), rootname, chainrec), "w")
for res in recinter :
    out_r.write("%s\n"%(res))

out_r.close()

out_l = open("%s/%s_%s_bs.txt"%(os.path.dirname(monol), rootname, chainlig), "w")
for res in liginter :
    out_l.write("%s\n"%(res))
out_l.close()

# writes res of bs in the bfactor
if bfact:

    if not generate : # must build the dict of rec and lig respectively
        dPDBlig = StructureTools.PDB_parser(monol)
        dPDBrec = StructureTools.PDB_parser(monor)

    write_interRes_inBfactor(dPDBrec, chainrec, recinter, dPDBlig, chainlig, liginter)
       
        
