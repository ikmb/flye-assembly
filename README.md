# flye-assembly
Assembles genomes from Pacbio reads using Flye

## This pipeline is currently designed to work on the Kiel MEDCluster (but could be adapted to other systems easily). 

## How to run

### Basic execution

`nextflow run ikmb/flye-assembly --bam /path/to/pacbio-movie.bam` --genome_size 1.2g

This requires Nextflow and Singularity. 

## Options

### `--bam`
A single movie file in BAM format (mutually exclusive with --samples). If you have multiple movies belonging to the same or different sequencing project/species, see the `--samples` option instead. 

### `--samples`
A CSV formatted sample sheet, using the following format:

`projectID;Movie;MovieIndex
MyGenome1;/path/to/movie.bam;/path/to/movie.bam.pbi
`

Add any number of movies to this list; all movies that share a project ID will be assembled together. Note that if you plan on assembling multiple genomes
using this approach, they need to share a genome size (at least roughly, i.e. within a factor of 2 or so). 

### `--genome_size`
The expected size of the genome(s) to be assembled. This option understands "g" and "m" to denote gigabases or megabases.

`--genome_size 1.3g` genome is expected to be 1.3 gigabases in size.

Note that the genome size only needs to be roughly accurate (within a factor of ~2).

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
