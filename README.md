# flye-assembly
Assembles genomes from Pacbio reads using Flye

## This pipeline is currently designed to work on the Kiel MEDCluster (but could be adapted to other systems easily). 

## How to run

### Basic execution

`nextflow run ikmb/flye-assembly --bam /path/to/pacbio-movie.bam`

## Options

### `--bam`
A movie file in BAM format (mutually exclusive with --samples)

### `--samples`
A CSV formatted sample sheet, using the following format:

`projectID;Movie;MovieIndex
MyGenome1;/path/to/movie.bam;/path/to/movie.bam.pbi
`

Add any number of movies to this list; all movies that share a project ID will be assembled together. Note that if you plan on assembling multiple genomes
using this approach, they need to share a genome size (at least roughly). 

### `--qc`
Switch on QC options - will generate information about the input reads and the resulting assembly

### `--qc_kat`
Run the KAT Kmer analysis 

### `--hifi`
Treat input movie as "hifi", and perform CCS generation prior to assembly

### `--reference`
This read data is from a previously sequenced species - use this reference genome for comparative analysis. 

* MiniMap read alignments and coverage analysis

### `--gff`
This read data is from a previously sequenced specues - this this reference gene annotation for comparative analysis.

* QUAST (requires --reference)
