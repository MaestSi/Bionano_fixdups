# Bionano_fixdups
**Bionano_fixdups** is a script for removing artificial duplications introduced by Bionano scaffolding pipeline. Starting from the identification of negative gaps annotated in the agp file, the script performs alignments between contigs at the 5' and at the 3' flanking regions of negative gaps in scaffolds, trims the overlaps, and produces a trimmed fasta file. The script is experimental, and its development was discontinued after the release of more refined tools as [BiSCoT](https://github.com/institut-de-genomique/biscot).

## Prerequisites

* Minimap2
* Samtools
* Jvarkit samextraclips
* R with BioStrings package

## Usage

**Bionano_fixdups.R**

Rscript ./Bionano_fixdups.R \<scaffolds.fasta\> \<file.agp\> \<contigs.fasta\>

*Note*: set the path to Minimap2, Samtools and Samextractlips executables before running the script.

Inputs:

*\<scaffolds.fasta\>: fasta file with scaffolds produced by Bionano hybrid scaffolding pipeline
*\<file.agp\>: agp file describing which contig has been included in each scaffold
*\<contigs.fasta\>: fasta file with contigs cut by Bionano hybrid scaffolding pipeline

Outputs:

* \<scaffolds\_neg\_gaps\_fixed.fasta\>: fasta file with overlaps between contigs trimmed
* logfile\_fix\_scaffolds.txt: logfile reporting operations performed on input scaffolds
