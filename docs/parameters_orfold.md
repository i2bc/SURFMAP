## ORFfold parameters


<b>Mandatory</b>

  ```-fna ```                 FASTA file containing the amino acid sequences to treat 



 ``` -options ```     indicates which properties are to be calculated. H for
estimating the fold potential with HCA, I for the estimation of the disorder
 propensity with IUPred and T for the aggregation propensity with Tango. Combinations
of letters are accepted if the user wants to calculate several properties at the
same time (-options HIT will estimate the three properties)(default: H).


<b>Optional</b>


  ```-h, --help ```           shows this help message and exits


 ```-gff```                  GFF annotation file. The ID (i.e. annotation) of the 
 sequences given in the input FASTA file sequence must be identical to the ID label 
 in the GFF file (column #3). ORFold generates as many GFF
files as studied properties (fold potential, disorder and/or aggregation 
 propensities), each containing for the sequences provided in the input FASTA file, 
 their corresponding property values (fold potential, disorder or aggregation 
 propensities). The values are stored in the column #9 of the output GFF files
 that can be subsequently uploaded on a genome viewer (see [here](./Run_orfold_advanced.md) 
 examples on the use of this option). 



```-keep``` ORFold uses IUPred and Tango for the prediction of the disorder 
and aggregation propensities. For storage reasons, by default, ORFold does not save the output 
files of these two methods. Nevertheless, the user can keep the Tango
output files through the option **-keep**. 
**-keep** T will save Tango output files in the TANGO directory.


```-N``` Working with large genomes can eventually generate quite big files. For
that reason, ORFold has also the option to generate a sample of the initial dataset 
and perform the calculation only on this specific subset of the population. This
option enables the user to indicate the number of randomly selected sequences to be treated by ORFold.

