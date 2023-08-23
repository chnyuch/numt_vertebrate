# numt_vertebrate
Scripts and software for NUMT identification and selection signals investigation.

## File descriptions
 - "numt_identification": scripts for NUMT, cNUMT identification
 - "numt_selection": scripts for dcNUMT clustering and selection signals identification

## Software requirement
Recommend to run on a multi-core workstation if many species are tested
### OS requirement
 - Linux: Ubuntu 20.04.6
### Python dependencies
```
os
sys
Bio
Bio.Seq
time
pandas
re
numpy
collections
ete3
glob
Bio.Phylo.PhyloXML
```
### Softwares requirement
- Genome download
  
  ncbi-genome-download

  https://github.com/kblin/ncbi-genome-download
  
- KaKs calculation
  
  KaKs_Calculator2.0

  https://sourceforge.net/projects/kakscalculator2/

- BLAST
  
  BLASTN v2.12.0+

  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/

- Selection signal investigation

  paml4.9j

  https://github.com/abacus-gene/paml

- Sequence extractions

  BEDTools v2.26.0

  https://github.com/arq5x/bedtools2

- Reverse sequence

  revseq in EMBOSS v6.6.0.0
  
  https://emboss.sourceforge.net/
 
- Sequence alignments

  MAFFT v7.310

  https://mafft.cbrc.jp/alignment/software/source.html

- Rate of evolution

  BayesTraits v.4.0.1

  https://www.evolution.reading.ac.uk/BayesTraitsV4.0.1/BayesTraitsV4.0.1.html

- cdNUMT and mito genes tree reconstructions

  iqtree v2.0.3

  https://github.com/iqtree/iqtree2

- Clustered NUMT tree recontructions

  FastTree v2.1.11

  http://www.microbesonline.org/fasttree/

-  Focal clade extraction

   Newick Utilities v1.7.0

   https://github.com/tjunier/newick_utils

## Script download
```
wget https://github.com/chnyuch/numt_vertebrate.git
```
### Instructions




  
