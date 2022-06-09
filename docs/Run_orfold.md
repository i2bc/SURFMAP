## Running ORFold:


### Inputs
Basically, ORFold requires only a FASTA file containing the amino acid sequences 
to treat (given with the **-fna** label). ORFold can handle several FASTA files at the same
time. In this case, it will treat them independently and will generate as many 
outputs as entered FASTA files.

 
 FASTA file example:
```{}
>aminoacid_sequence_1
AGNVCFGGRTYMPDFDGMSCVNWQERT
>aminoacid_sequence_2
MPDFMPCNVSDRTEEEPMSPARTYDFGHKLCVSDFTPMLKKPERT
```


### How to estimate the fold potential and/or disorder and aggregation propensities
By default, ORFold only estimates the fold potential of the input sequences. 
The disorder and aggregation propensities can be however calculated as well.
The user can specify which calculation methods are to be launched with 
the **-options** argument. 

Each method used by ORFold is referred by its initial: 
<pre>
   HCA     : H
   IUPred  : I
   TANGO   : T 
</pre>

The user must specify the combination of methods he wants to apply
on the input sequences giving their initials with the **-options** argument
without any space: ```-options HIT``` for running the 3 programs or ```-options HT``` 
if the user wants to run only HCA and Tango for example. 
The order of the letters has no importance, 
```-options HIT``` and ```-options THI``` will lead to the same result.


### Basic run
The following instruction estimates the fold potential, and the disorder and aggregation propensities of
all amino acid sequences contained in the input fASTA file:

```{bash}
orfold -fna sequences.fasta -options HIT
```



The user has to notice that **IUPred** and **Tango** provide additional information
to HCA but will slow down considerably ORFold for large datasets. 
The next instruction only calculate the fold potential with HCA:
```{bash}
orfold -fna sequences.fasta -options H
```

### Output:
ORFold produces a table (fasta_rootname.tab) that contains for each input sequence,
the computed values (separated by tabulations) according to the user request (fold potential, and/or disorder and/or 
aggregation propensities).





 Output file example with ```-options HIT``` (fold potential, disorder and aggregation
 propensities estimated from HCA, IUPred and Tango, see [here](./How_it_works_orfold.md) for more details):

```{}
Seq_ID                  HCA     Disord  Aggreg
aminoacid_sequence_1	1.340	0.000	0.230	
aminoacid_sequence_2	-0.230	0.120	0.012	
```

Output file example with ```-options H``` (fold potential estimated with HCA, see [here](./How_it_works_orfold.md) for more details):
```{}
Seq_ID                  HCA     Disord  Aggreg
aminoacid_sequence_1    1.340   nan     nan
aminoacid_sequence_2    -0.230  nan     nan
```
