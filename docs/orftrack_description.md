![LOGO_ORFtrack](./img/icons/Logo_ORFtrack.png){ width=30% }
## Aims and general description

ORFtrack scans a given genome in the six frames, and searches for 
all possible ORFs longer than a given size (default: 60 nucleotides -
STOP codons excluded). It annotates them according to a set of genomic features (e.g. noncoding intergenic,
coding, noncoding and overlapping with a specific genomic feature - see
the [ORF annotation section](./orftrack_annotation.md) for more details). 
ORFtrack takes as inputs a FASTA file containing the nucleotide
sequences of all chromosomes or contigs and their corresponding 
annotations in a GFF file. The program returns a new GFF file that contains all
identified ORFs. In addition, the amino acid sequences of 
all annotated ORFs or specific subsets of ORFs (i.e. only noncoding intergenic ORFs for example)
can be extracted and written in a FASTA file with ORFget, a tool 
provided with ORFtrack. 
