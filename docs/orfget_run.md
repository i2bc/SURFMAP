## Extraction and writing of ORF sequences with ORFget

ORFget is a tool provided with ORFtrack that allows the user to extract
the protein and/or nulceotide sequences of specific subsets of ORFs 
according to their annotation categories
(see [here](./orftrack_annotation.md) 
for a description of all ORF categories). ORFget deals with annotation 
patterns, thereby allowing different levels of annotation in a 
very easy fashion.

ORFget has two principal options:

* ```-features_include```: list of motifs that will be used to define the 
  ORFs that will be included in the FASTA 
  output. The sequences whose annotations include these patterns will 
  be retained in the output FASTA file 
* ```-features_exclude```: list of motifs that will be used to define the 
  ORFs that will be excluded in the FASTA 
  output. The sequences whose annotations include these patterns 
  will not be written in the output FASTA file
  
The searched patterns can be specific (for a finer selection) or more general.<br><br>

<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
For example, the motif <b>"nc"</b> appears in the features: <b>nc</b>_intergenic, <b>nc</b>_ovp_same-mRNA, <b>nc</b>_ovp_opp-mRNA and <b>nc</b>_ovp_same-tRNA.
As a result, the option <code>-feature_include nc</code> will keep all the four
feature categories. 
<br><br>
The option <code>-feature_include nc_ovp</code> will keep:
<ul>	
 <li><b>nc_ovp</b>_same-mRNA</li>
 <li><b>nc_ovp</b>_opp-mRNA</li>
 <li><b>nc_ovp</b>_same-tRNA</li>
</ul>

The option <code>-feature_include nc_ovp_same</code> will keep:
<ul>
 <li><b>nc_ovp_same</b>-mRNA</li>
 <li><b>nc_ovp_same</b>-tRNA</li>
</ul>

The option <code>-feature_include mRNA</code> will keep: 
<ul>
 <li>nc_ovp_same-<b>mRNA</b></li>
 <li>nc_ovp_opp-<b>mRNA</b></li>
</ul>

The option <code>-feature_exclude opp</code> will eliminate the nc_ovp_<b>opp</b>-mRNA and will keep:
<ul>
 <li>nc_intergenic</li>
 <li>nc_ovp_same-mRNA</li>
 <li>nc_ovp_same-tRNA</li>
etc... 
</p>
</div>
Here are presented some examples of selection of ORFs with ORFget.


### Extraction of the sequences of all the ORFs of a GFF file

The following command writes the amino acid sequences of all ORFs 
annotated in the input GFF file.


``` python
orfget -fna genome.fasta -gff mapping_orf_genome.gff
```
ORFget generates a FASTA file containing all the corresponding amino acid
sequences. 



### Extraction of the sequences of all noncoding ORFs identified with ORFtrack

The following commands, each enable the user to write the 
amino acid sequences of all noncoding 
ORFs no matter their status (i.e. intergenic or overlapping)
(see [here](./orftrack_annotation.md) for a description of all ORF categories).

``` bash
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_include nc
```
or 
``` bash
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_include nc_intergenic nc_ovp
```
or
``` bash
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_exclude c_CDS
```

### Extraction of the sequences of a specific subset of ORFs according to their annotation

The following instruction writes the amino acid sequences of the ORFs
which overlap with CDS on the same, or on the opposite strand.

``` bash
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_include nc_ovp_same-CDS nc_ovp_opp-CDS
```


Notice that using the argument "features_exclude" assumes that the selection 
operates on all genomic features except those that are excluded. 
Consequently, if the user wants to select all noncoding sequences
except those overlapping CDS, mRNAs, tRNAs, and rRNAs, he must 
exclude the coding ORFs (c_CDS) as well. Otherwise, they will be
kept.


``` bash
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_exclude c_CDS nc_same_ovp-tRNA nc_same_ovp-rRNA nc_opp_ovp-mRNA nc_opp_ovp-tRNA nc_opp_ovp-rRNA nc_opp_ovp-mRNA  
```

### Extraction of the sequences of a random subset of ORFs 

Sometimes, for computational time or storage reasons, the user does 
not want to deal with all the ORFs of a specific category. ORFget
can provide the user with a subset of N (to be defined by the user)
randomly selected ORFs from a specific ORF category. The last instruction
writes the sequences of 10000 randomly selected noncoding 
intergenic ORFs.


``` python
orfget -fna genome.fasta -gff mapping_orf_genome.gff -features_include nc_intergenic -n 10000
```

### Reconstruction of protein sequences
In addition, ORFget enables the reconstruction of all protein 
sequences of a genome (i.e. all isoforms) according to their 
definition in the original GFF file. The following instruction
writes all the resulting sequences in a FASTA file.


``` python
orfget -fna genome.fasta -gff genome.gff -features_include CDS
```

### Writing amino acid or nucleotide sequences
By default, ORFget will generate the amino acid sequences of the 
desired ORFs in a FASTA file 
with the extension **.pfasta**. If the user wishes to generate the nucleotide
or even both nucleotide and amino acids sequences, he must use the 
option
```-type nucl``` and ```-type both```, respectively.

