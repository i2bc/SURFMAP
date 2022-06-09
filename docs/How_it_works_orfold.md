# How works ORFold?

ORFold is a tool developed in python3 which aims at 
characterizing the fold potential of a set of amino acid sequences with no 
knowledge of their 3D structures nor evolutionary information (orphan sequences can be treated). 
The fold potential of a sequence is 
calculated with the HCA 
method. Also, ORFold can estimate the disorder, 
and the aggregation propensities of the input sequences with IUPred
and Tango respectively.    


## Hydrophobic Clusters Analysis (HCA)
**HCA**[1] aims as delineating in an amino acid sequence, regions enriched in 
strong hydrophobic residues (HCA clusters) and regions 
of at least four consecutive non-hydrophobic residues (HCA linkers). 
The patterns of hydrophobic residues can be associated with specific regular 
secondary structures, and the distribution of the HCA clusters and linkers in a protein 
sequence can be used to estimate through the HCA score, its ability to fold (completely or partially). 

This score ranges from -10 to +10 with low HCA scores indicating
sequences depleted in hydrophobic clusters and expected to be disordered in solution, 
while high HCA scores reflect sequences enriched in hydrophobic clusters 
and expected to generate aggregates in solution, though some of them could
fold in lipidic environments. 
Foldable sequences are known to display
an equilibrium between hydrophobic and hydrophilic residues (average of 33% 
of hydrophobic residues in globular proteins). 
They are mostly associated with intermediate HCA score values.
 

The HCA score is calculated using the freely available 
software **pyHCA** which can be downloaded and installed 
following the instructions of its developers: <https://github.com/T-B-F/pyHCA>


## Tango
**Tango**[5][6][7] is a method which aims at predicting aggregation nucleating regions
in protein sequences. 
If specified by the user, ORFold can calculate and add the aggregation propensity 
of a sequence in the output. 
Tango is not freely available software, and the user of ORFold should 
first contact the Tango developers to have access to the source code: <http://tango.crg.es>

For the aggregation propensity estimation, according to the protocol
proposed by Linding et al.[6], a sequence segment is considered as aggregation prone
if it is composed of at least five consecutive residues predicted as 
populating a b-aggregated conformation with a percentage occupancy greater than 5%. 
Then, the aggregation propensity of each sequence is defined as the 
fraction of residues predicted in aggregation prone segments. 

## IUPred
**IUPred2A**[2][3][4] is one of the best methods for the prediction of 
Intrinsically Disordered Proteins (IDPs) and can be used as a 
complement to the HCA score prediction. 
If specified by the user, ORFold can calculate and add the disorder propensity 
of a sequence in the output. 
IUPred is not freely available, and the user of 
ORFold should first contact the IUPred developers to 
have access to the source code : <https://iupred2a.elte.hu>

For the disorder propensity estimation, in order to be consistent with
the estimation of the aggregation propensity, ORFold searches for 
regions on the protein sequence that present at least five consecutive 
residues with a disorder probability higher than 0.5. 
The disorder propensity of each sequence is defined as the fraction 
of residues predicted as located in a highly disordered segment.    



<br><br><br>
#### References

1. Bitard-Feildel, T. & Callebaut, I. HCAtk and pyHCA: A Toolkit and Python API for the Hydrophobic Cluster Analysis of Protein Sequences. bioRxiv 249995 (2018).
2. Dosztanyi, Z., Csizmok, V., Tompa, P. & Simon, I. The pairwise energy content estimated from amino acid composition discriminates between folded and intrinsically unstructured proteins. Journal of molecular biology 347, 827–839 (2005).
3. Dosztányi, Z. Prediction of protein disorder based on IUPred. Protein Science 27, 331– 340 (2018).
4. Mészáros, B., Erdős, G. & Dosztányi, Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic acids research 46, W329–W337 (2018).
5. Fernandez-Escamilla, A.-M., Rousseau, F., Schymkowitz, J. & Serrano, L. Prediction of sequence-dependent and mutational effects on the aggregation of peptides and proteins. Nature biotechnology 22, 1302–1306 (2004).
6. Linding, R., Schymkowitz, J., Rousseau, F., Diella, F. & Serrano, L. A comparative study of the relationship between protein structure and β-aggregation in globular and intrinsically disordered proteins. Journal of molecular biology 342, 345–353 (2004).
7. Rousseau, F., Schymkowitz, J. & Serrano, L. Protein aggregation and amyloidosis: confusion of the kinds? Current opinion in structural biology 16, 118–126 (2006).
