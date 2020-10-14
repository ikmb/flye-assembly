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

Mandatory arguments (one of):

--samples		A sample sheet with information on project and location of associated movies
OR
--bam			A movie file in BAM format (must be accompanied by a .pbi index file)

Options:
--hifi			Wether the reads should be treated as hifi (i.e. ultra-high subread coverage)
--genome_size		Expected genome size 
--qc			Wether to run quality control steps
--qc_kat		Enable KAT kmer analysis (may crash...)
--reference		A reference genome to compare against (also requires --gff)
--gff			A reference annotation (also requires --reference)
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
log.info "=================================================="
if (params.samples) {
	log.info "SampleSheet:			${params.samples}"
}
if (params.bam) {
	log.info "Movie file:		${params.bam}"
}
log.info "IsHifi:			${params.hifi}"
log.info "Genome size:		${params.genome_size}"
log.info "Run QC:			${params.qc}"
log.info "Run KMER Analaysis:	${params.qc_kat}"
if (workflow.containerEngine) {
        log.info "Container engine:     	${workflow.containerEngine}"
}
log.info "================================================="
log.info "Starting at:          $workflow.start"

if (params.reference && params.gff) {
	log.info "Will run QUAST for quality control!"
}
		
// Enables splitting of CCS read generation into 10 parallel processes
def chunks = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]

if (params.reference) {
	Ref = Channel.fromPath(file(params.reference))
} else {
	Ref = Channel.empty()
}

def returnFile(it) {
    // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "Missing file in TSV file: ${inputFile}, see --help for more information"
    return inputFile
}

if (params.bam) {
	bamFile = Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "Could not find an input bam file" }
		.map { b -> [ file(b).getBaseName(), file(b), file("${b}.pbi") ] }
} else if (params.samples) {
	// Read sample file
	bamFile = Channel.from(file(params.samples))
	       	.splitCsv(sep: ';', header: true)
		.map { row -> 
			def project_name = row.ProjectID
			def movie = returnFile( row.Movie )
			def movie_index = returnFile( row.MovieIndex )
			[ project_name, movie, movie_index ]
		}

} else {
        exit 1, "Must provide a sample sheet or a movie in BAM format (see documentation for details)"
}

if (!params.genome_size) {
	exit 1, "Must provide an approximate genome size (--genome_size)"
}


/*
Start the pipeline here
*/

if (params.hifi) {

	process BamToCCS {

		publishDir "${params.outdir}/${sample}/CCS", mode: 'copy',
			saveAs: { filename ->
				if( filename.indexOf("report.txt") > 0 ) filename
				else null
			}

		scratch true 

		label 'pbccs'

		input:
		set val(sample), file(bam), file(bam_index) from bamFile
		each chunk from chunks

		output:
		set val(sample),file(reads) into ReadChunks
		file(report)

		script:
		reads = bam.getBaseName() + "." + chunk + ".ccs.bam"
		report = bam.getBaseName() + "." + chunk + ".ccs.report.txt"

		"""
			ccs $bam $reads --chunk $chunk/10 -j ${task.cpus}
			mv ccs_report.txt $report
		"""

	}

	ReadChunksGrouped = ReadChunks.groupTuple()

	process CcsMerge {

		label 'pbbam'

		publishDir "${params.outdir}/${sample}/CCS", mode: 'copy'

		input:
		set val(sample),file(read_chunks) from ReadChunksGrouped

		output:
		set val(sample),file(bam),file(pbi) into mergedReads

		script:
		bam = sample + ".reads.ccs.bam"
		pbi = bam + ".pbi"

		"""
			pbmerge -o $bam $read_chunks
			pbindex $bam
		"""
	}

} else {

	mergedReads = bamFile

}

process BamToFastq {

	label 'bam2fastx'

	publishDir "${params.outdir}/${sample}/fastq", mode: 'copy'

	scratch true

        input:
        set val(sample),file(bam),file(pbi) from mergedReads

        output:
        set val(sample),file(reads) into (Reads, ReadsNanoqc, ReadsKat, ReadsNanoplot, ReadsAlign)

        script:
        reads = bam.getBaseName() + ".fastq.gz"

        """
        	bam2fastq -o ${bam.getBaseName()} $bam
        """

}

grouped_movies = Reads.groupTuple()

process Flye {

	publishDir "${params.outdir}/${sample}/assembly", mode: 'copy'

	label 'flye'

	input:
	set val(sample),file(reads) from grouped_movies

	output:
	set val(sample),file(assembly) into ( Assembly, AssemblyBusco )
	file(assembly_renamed) into AssemblyQuast
	file(assembly_info)

	script:
	folder_name = "flye_assembly"
	run_options = ""
	if (params.hifi) {
		folder_name = "flye_assembly_hifi"
	}
	assembly = folder_name + "/assembly.fasta"
	assembly_info = folder_name + "/assembly_info.txt"
	assembly_renamed = sample + ".assembly.fasta"

	def options
	if (params.hifi) {
		options = "--pacbio-hifi"
	} else {
		options = "--asm-coverage 40 --pacbio-raw"
	}

	"""
		flye $options $reads -i 2 --genome-size ${params.genome_size} --threads ${task.cpus} --out-dir $folder_name
		cp $assembly $assembly_renamed
	"""

}

process nanoQC {

	label 'nanoqc'

	publishDir "${params.outdir}/${sample}/QC", mode: 'copy'

	when:
	params.qc

	input:
	set val(sample),file(fastq) from ReadsNanoqc

	output:
	file(qc_plot)

	script:
	qc_plot = "nanoqc/nanoQC.html"

	"""
		nanoQC -l 100 -o nanoqc $fastq
	"""
}

process nanoPlot {

	label 'nanoplot'

        publishDir "${params.outdir}/${sample}/QC", mode: 'copy'

	when:
	params.qc

	input:
	set val(sample),file(fastq) from ReadsNanoplot

	output:
	file("nanoplot")

	script:
	"""
		NanoPlot -t ${task.cpus} -o nanoplot --fastq $fastq 
	"""	
}

process QCN50 {

	publishDir "${params.outdir}/${sample}/QC", mode: 'copy'

	label 'gaas'

	when:
	params.qc

	input:
	set val(sample),file(assembly) from Assembly

	output:
	file(assembly_report)

	script:
	assembly_report = "stats/" + assembly.getBaseName() + "_stat.txt"
	"""
		gaas_fasta_statistics.pl -f $assembly -o stats
	"""

}

process KatGCP {

        publishDir "${params.outdir}/${sample}/QC", mode: 'copy'

	label 'kat'

	when:
	params.qc_kat

	input:
	set val(sample),file(reads) from ReadsKat

	output:
	file(kat_dist) into KatDist

	script:
	kat_dist = "kat-gcp.mx.png"

	"""
		kat gcp -t ${task.cpus} $reads
	"""

}

process Quast {

	publishDir "${params.outdir}/Quast", mode: 'copy'

	label 'quast' 

	when:
	params.reference && params.gff

	input:
	file(assemblies) from AssemblyQuast.collect()

	output:
	file(quast_dir)

	script:
	quast_dir = "quast"
	"""
		quast -o $quast_dir -r ${params.reference} --features ${params.gff} --threads ${task.cpus} $assemblies
	"""

}

process makeRef {

	label 'map'

	when:
	params.reference

	input:
	file(fasta) from Ref

	output:
	set file(idx),file(fasta) into RefMap

	script:
	idx = fasta.getBaseName() + ".idx"

	"""
		minimap2 -d $idx $fasta
	"""
}

process Minimap {

	label 'map'

	publishDir "${params.outdir}/${sample}/Mapping", mode: 'copy'

	when:
	params.reference 

	input:
	set val(sample),file(reads) from ReadsAlign
	set file(idx),file(fasta) from RefMap.collect()

	output:
	set val(sample),file(bam),file(bai)

	script:
	bam = reads.getBaseName() + ".bam"
	bai = bam + ".bai"

	"""
		minimap2 -ax map-pb $idx $reads | samtools sort -m 8G - | samtools view -bh -o $bam -
		samtools index $bam
	"""

}

