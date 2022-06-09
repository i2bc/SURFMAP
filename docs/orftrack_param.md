## ORFtrack parameters


<b>Mandatory</b>

  ```-fna ```                 nucleotide fasta file of the genome whose ORFs are to annotate

  ```-gff```                  GFF annotation file


<b>Optional</b>


  ```-h, --help```            shows the help and exits

 ```` --show-types ````       prints all genomic features annotated in the input GFF file
 
  ```--show-chrs ```          prints all the seqID (usually corresponding to chromosome or contig ID) present in the input GFF file


 ```` -chr````                  list of seqID to be treated by ORFtrack (i.e. column #1 
                        of the GFF file - generally chromosome or contig ID). 
                        ID must be separated by a space: -chr NC_001148.4 NC_001139.3   
In this case, ORFtrack will not treat the entire genome, but
the following seqID NC_001148.4 NC_001139.3. We recommend using this option 
when dealing with large genomes in order to distribute the calculations on 
several CPUs (one per (subset of) seqID).



  ```-orf_len ```  Minimal number of nucleotides between two 
  consecutive STOP codons to define an ORF (default: 60 nucleotides) 
(see the [ORF definition section](./orftrack_orfdef.md) for more details).

 ```` -co_ovp````  Minimal fraction of the ORF length that overlaps a genomic feature
 to annotate the ORF as overlapping it 
 (default: 70%). (see the [Overlap section](./orftrack_overlap.md) for more details).


  ```-out```           Output directory



  `````-types_only````` Genomic feature(s) considered for the annotation step ('CDS' is
                        included by default). If there are several genomic features
  to be considered, they must be given separated by a space. Noncoding ORFs are annotated as
intergenic (when they do not overlap any feature) or overlapping (when 
  they overlap a given genomic feature). The "overlapping" status 
directly derives from the genomic features considered for the annotation
step. For example, if the user specifies with the "types_only" option, the features
"tRNA" and "rRNA" (CDS are included by default), all ORFs that overlap another
  genomic feature (i.e. different from tRNA, rRNA, or CDS) will be annotated as 
  noncoding intergenic ORFs. If no genomic feature is indicated, noncoding ORFs that overlap 
  any genomic feature annotated in the original GFF file will be annotated 
 as noncoding ORF overlapping the corresponding genomic feature. Nevertheless,
the resulting ORFs can be treated specifically (grouped according to their annotation
  for example) with ORFget afterwards
  (see the [ORFget section](./orfget_run.md) for more details).


  ```-types_except```   Genomic feature(s) not to be considered for the annotation step 
                        ('gene' and 'exon' not considered by default). If there are several genomic features
  to be considered, they must be given separated by a space. Noncoding ORFs are annotated as
intergenic (when they do not overlap any feature) or overlapping (when 
  they overlap a given genomic feature). The "overlapping" status 
directly derives from the genomic features considered for the annotation
step. For example, when specifying with the "types_except" option, 
  the features "telomer" and "centromer" ('gene' and 'exon' not considered by 
  default), if an ORF overlaps a "telomer", "centromer", "gene" and/or "exon", 
the overlap will not be considered and the ORF will be annotated as intergenic or 
  overalapping another genomic feature if there is another genomic feature 
in the same region (e.g. an ORF that overlaps at the same time a gene and a 
  tRNA will be annotated as overlapping a tRNA since the gene is not considered
  for the ORF annotation). If no genomic feature is indicated, noncoding ORFs that overlap 
  any genomic feature annotated in the original GFF file will be annotated 
 as noncoding ORF overlapping the corresponding genomic feature. Nevertheless,
the resulting ORFs can be treated specifically (grouped according to their annotation
  for example) with ORFget afterwards
  (see the [ORFget section](./orfget_run.md) for more details).

 