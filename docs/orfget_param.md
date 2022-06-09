## ORFget parameters


<b>Mandatory</b>

  ```-fna ```                 nucleotide fasta file of the genome whose ORFs are to extract


 ```-gff```                  GFF annotation file


  ````-o````                 root name of the output fasta

<b>Optional</b>


  ```-h, --help ```           shows this help message and exits

 
  ```-features_include```  genomic features or annotation patterns 
                           to be searched for the ORF selection.
ORFs whose annotations match the entered patterns or genomic features 
are written in the output fasta file. (default: all genomic features present in 
  the input GFF, consequently, all annotated ORFs will be written in the output
GFF file)(see the [ORFget section](./orfget_run.md) for more details).


  ```-features_exclude```  genomic features or annotation patterns 
                           to be searched for the ORF selection.
ORFs whose annotations match the entered patterns or genomic features 
are not written in the output fasta file. (default: None, consequently, 
  all annotated ORFs will be written in the output GFF file)(see the
  [ORFget section](./orfget_run.md) for more details).



  ```-chr_exclude```  to exclude a (or a subset of) seqID(s) for the ORF
extraction and writing - the ORFs belonging to this seqID (usually a chromosome
  or a contig) will not be written in the output FASTA file no matter their
  annotation category (default: None). 


  ```-N```       number of ORFs randomly selected for the writing according
to their genomic feature or pattern.


