# flye-assembly
Assembles genomes from Pacbio reads using Flye

## This pipeline is currently designed to work on the Kiel MEDCluster (but could be adapted to other systems easily). 

## How to run

### Basic execution

`nextflow run ikmb/flye-assembly --bam /path/to/pacbio-movie.bam`

## Options

### `--qc`
Switch on QC options - will generate information about the input reads and the resulting assembly

### `--qc_kat`
Run the KAT Kmer analysis 

### `--hifi`
Treat input movie as "hifi", and perform CCS generation prior to assembly

### `--busco_lineage`
Switches on BUSCO to gauge gene space coverage and specifies which reference dataset will be used. Available options:

* actinopterygii_odb10 
* carnivora_odb10
* cyprinodontiformes_odb10
* eutheria_odb10
* glires_odb10
* mammalia_odb10
* mollusca_odb10
* passeriformes_odb10
* sauropsida_odb10
* vertebrata_odb10
* aves_odb10
* cetartiodactyla_odb10
* euarchontoglires_odb10
* laurasiatheria_odb10
* metazoa_odb10
* nematoda_odb10
* primates_odb10
* tetrapoda_odb10

### `--reference`
This read data is from a previously sequenced species - use this reference genome for comparative analysis. 

* MiniMap read alignments and coverage analysis

### `--gff`
This read data is from a previously sequenced specues - this this reference gene annotation for comparative analysis.

* QUAST (requires --reference)
