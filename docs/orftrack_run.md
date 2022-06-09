## Annotation of ORFs with ORFtrack

### Running ORFtrack on a complete genome

The following instruction annotates all the possible ORFs of the
input genome.


``` bash
orftrack -fna genome.fasta -gff genome.gff
```
Depending on the size of the genome, ORFtrack takes a few minutes to
several hours to annotate all the ORFs of a genome. It returns 
a new GFF file containing the annotation of all the identified ORFs 
(including coding and noncoding ORFs). 


### Running ORFtrack on a single chromosome or a subset of chromosomes

ORFtrack can be launched on a single seqID (usually chromosome or contig indicated in the first 
column of the inout GFF)(e.g. chromosome seqID: XXX) 
with the following instruction:


``` bash
orftrack -fna genome.fasta -gff genome.gff  -chr chr_ID_XXXX
```
This can be very useful if the user wants to run ORFtrack on several 
CPUs. Also, it can be launched on a subset of seqIDs as follows:


``` bash
orftrack -fna genome.fasta -gff genome.gff  -chr seqID1 seqID2 seqIDx
```



