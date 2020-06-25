#!/usr/bin/env nextflow

/*

||  ||
Flye Assembly
||  ||

*/

def helpMessage() {
  log.info"""
======================================
>> Flye assemmbly from PacBio reads <<
======================================

Usage:

A typical command to run this pipeline would be:

nextflow run ikmb/flye-assembly --bam movie.bam --genome_size 1.3g --qc

Mandatory arguments:

--bam			A PacBio movie file in BAM format

Options:
--hifi			Wether the reads should be treated as hifi (i.e. ultra-high subread coverage)
--genome_size		Expected genome size 
--qc			Wether to run quality control steps
--email                 Provide an Email to which reports are send.
--skip_multiqc          Do not generate a QC report
""".stripIndent()
}

// Show help message
if (params.help){
        helpMessage()
        exit 0
}

log.info "=================================================="
log.info "Flye assembly pipeline v${workflow.manifest.version}"
log.info "Nextflow Version:     $workflow.nextflow.version"
log.info "Command Line:         $workflow.commandLine"
log.info "Authors:              M. Höppner"
if (workflow.containerEngine) {
        log.info "Container engine:     ${workflow.containerEngine}"
}
log.info "================================================="
log.info "Starting at:          $workflow.start"


def chunks = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]

if (!params.bam) {
	exit 1, "Must provide a BAM movie file (--bam)"
} else {
	bamFile = Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "Could not find an input bam file" }
		.map { b -> [ file(b), file("${b}.pbi") ] }
}

if (!params.genome_size) {
	exit 1, "Must provide an approximate genome size (--genome_size)"
}

/*
Start the pipeline here
*/


if (params.hifi) {

	process BamToCCS {

		label 'pbccs'

		input:
		set file(bam),file(bam_index) from bamFile
		each chunk from chunks


		output:
		file(reads) into ReadChunks

		script:
		reads = bam.getBaseName() + ".ccs.bam"

		"""
			ccs $bam $reads --chunk $chunk/10 -j ${task.cpus}
		"""

	}

	process CcsMerge {

		label 'samtools'

		input:
		file(read_chunks) from ReadChunks.collect()

		output:
		file(bam) into mergedReads

		script:
		reads = "reads.ccs.bam"

		"""
			samtools merge -@ ${task.cpus} $bam $read_chunks
		"""
	}

	process CcsToFastq {

		label 'bam2fastx'

		input:
		file(bam) from mergedReads

		output:
		file(reads) into Reads

		script:
		reads = bam.getBaseName() + ".fastq.gz"

		"""
			
		"""

	}

} else {

	process BamToFastq {

		label 'bam2fastx'

		input:
		file(bam) from bamFile

		output:
		file(reads) into Reads

		script:
		reads = bam.getBaseName() + ".fastq.gz"

		"""
			bam2fastq -o ${bam.getBaseName()} $bam
		"""

	}

}

process Flye {

	label 'flye'

	input:
	file(reads) from Reads

	output:
	file(assembly) into Assembly
	file(assembly_info)

	script:
	assembly = "flye_assembly/assembly.fasta"
	assembly_info = "flye_assembly/assembly_info.txt"

	def options
	if (params.hifi) {
		options = "--pacbio-hifi"
	} else {
		options = "--pacbio-raw"
	}

	"""
		flye $options $reads --genome-size ${params.genome_size} --threads ${task.cpus} --out-dir flye_assembly
	"""

}
