# Bionano_fixdups
**Bionano_fixdups** is a script for removing artificial duplications introduced by Bionano scaffolding pipeline. The project was 

## Prerequisites

* Minimap2
* Samtools
* Jvarkit samextraclips
* R with BioStrings package

## Usage

**Bionano_fixdups.R**

Rscript ./Bionano_fixdups.R \<scaffolds\> \<agp\> \<contigs\>
  
Inputs:
\<scaffolds\>: fasta file with scaffolds produced by Bionano hybrid scaffolding pipeline
\<agp\>: agp file describing which contig has been included in each scaffold
\<contigs\>: fasta file with contigs cut by Bionano hybrid scaffolding pipeline

Outputs:

**Note**: before running, edit the script to set variables SAMEXTRACTCLIPS, MINIMAP2, SAMTOOLS and SEQTK to the corresponding executables


Note: modify **config_MetaBlast.sh** before running; 

Inputs:
