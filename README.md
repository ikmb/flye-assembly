# flye-assembly
Assembles genomes from Pacbio reads using Flye

## This pipeline is currently designed to work on the Kiel MEDCluster (but could be adapted to other systems easily). 

## How to run

### Basic execution

`nextflow run ikmb/flye-assembly --bam /path/to/pacbio-movie.bam --genome_size 1.2g`

This requires Nextflow and Singularity. 

## Options

### `--bam`
A single movie file in BAM format (mutually exclusive with --samples). If you have multiple movies belonging to the same or different sequencing project/species, see the `--samples` option instead. 

### `--samples`
A CSV formatted sample sheet, using the following format:

```
projectID;Movie;MovieIndex

MyGenome1;/path/to/movie.bam;/path/to/movie.bam.pbi
```

Add any number of movies to this list; all movies that share a project ID will be assembled together. Note that if you plan on assembling multiple genomes
using this approach, they need to share a genome size (at least roughly, i.e. within a factor of 2 or so). 

### `--flye`
Run the Flye assembler

### `--hifiasm` 
Run the HifiASM assembler (requires --hifi)

### `--ipa`
Run the IP assembler (requires --hifi)

### `--large`
Boolean flag to enable assembly of larger genomes. By default, flye will run on a normal "fat" compute node. With this option it is possible to divert it to a dedicated "high-memory" node that may have specs different
from the standard compute nodes. On Medcluster, this will enable the use of 1.5TB Ram and 30 cores on the "assembly node". 

### `--genome_size`
The expected size of the genome(s) to be assembled. Specifying a genome size is optional, but if you do, Flye can downsample the reads during assembly to reduce compute requirements. 
This option understands "g" and "m" to denote gigabases or megabases. This information does NOT currently inform the aforementioned option `--large`. 

`--genome_size 1.3g` genome is expected to be 1.3 gigabases in size.

Note that the genome size only needs to be roughly accurate (within a factor of ~2). 

### `--qc`
Switch on QC options - will generate information about the input reads and the resulting assembly

### `--qc_kat`
Run the KAT Kmer analysis 

### `--trimming`
Remove 15bp from both ends of each read. Do this, if you suspect that reads retain adapter/barcode fragments. 

### `--hifi`
Treat input movie as "hifi", and perform CCS generation prior to assembly

### `--reference`
This read data is from a previously sequenced species - use this reference genome for comparative analysis. 

* MiniMap read alignments and coverage analysis

### `--gff`
This read data is from a previously sequenced specues - this this reference gene annotation for comparative analysis.

* QUAST (requires --reference)

### `-profile`
The execution profile to use. Defaults to "standard", which corresponds to the Medcluster. 

* standard  - CAU Medcluster
* assembly - IKMB assembly machine


