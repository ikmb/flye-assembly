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
--qc_kat		Enable KAT kmer analysis (may crash...)
--reference		A reference genome to compare against (also requires --gff)
--gff			A reference annotation (also requires --reference)
--busco			Run BUSCO assembly
--busco_lineage		Busco profile to use (e.g. aves_odb10)
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
log.info "Movie:			${params.bam}"
log.info "IsHifi:			${params.hifi}"
log.info "Genome size:		${params.genome_size}"
log.info "Run QC:			${params.qc}"
log.info "Run KMER Analaysis:	${params.qc_kat}"
if (params.busco_lineage) {
	log.info "BUSCO databases:	${params.busco_database_dir}"
	log.info "BUSCO lineage:		${params.busco_lineage}"
} else {
	log.info "Run BUSCO:		${params.busco}"
}
if (workflow.containerEngine) {
        log.info "Container engine:     	${workflow.containerEngine}"
}
log.info "================================================="
log.info "Starting at:          $workflow.start"

if (params.reference && params.gff) {
	log.info "Will run QUAST for quality control!"
}

if (params.busco_lineage) {
	if (!params.busco_database_dir) {
		exit 1, "Cannot run Busco without a database directory (--busco_database_dir)"
	}
	lineage_dir = file("${params.busco_database_dir}/${params.busco_lineage}")
	if (!lineage_dir.exists()) {
		exit 1, "Did not find the specified database dir / lineage on this machine"
	}
	println "Will run BUSCO - make sure you have Miniconda loaded or this will fail...!"
}
		
// Enables splitting of CCS read generation into 10 parallel processes
def chunks = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]

if (params.reference) {
	Ref = Channel.fromPath(file(params.reference))
} else {
	Ref = Channel.empty()
}

if (!params.bam) {
	exit 1, "Must provide a BAM movie file (--bam)"
} else {
	bamFile = Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "Could not find an input bam file" }
		.map { b -> [ file(b).getBaseName(), file(b), file("${b}.pbi") ] }
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
		set val(sample),file(bam),file(bam_index) from bamFile
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

		publishDir "${params.outdir}/${sample}/CCS", mode: 'copy',

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

	mergedReads = bamFile.groupTuple()

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

process Flye {

	publishDir "${params.outdir}/${sample}/assembly", mode: 'copy'

	label 'flye'

	input:
	set val(sample),file(reads) from Reads

	output:
	set val(sample),file(assembly) into ( Assembly, AssemblyBusco )
	file(assembly_renamed) into AssemblyQuast
	file(assembly_info)

	script:
	folder_name = "flye_assembly"
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
		options = "--pacbio-raw"
	}

	"""
		flye $options $reads -i 2 --genome-size ${params.genome_size} --threads ${task.cpus} --out-dir $folder_name
		cp $assembly $assembly_renamed
	"""

}

process Busco {

	publishDir "${params.outdir}/${sample}/Busco", mode: 'copy'

	label 'busco'

	when:
	params.busco_lineage

	input:
	set val(sample),file(assembly) from AssemblyBusco
	
	output:
	file(busco_report)

	script:
	busco_report = "busco/short_summary.specific." + params.busco_lineage + ".busco.txt"

	"""
		busco -i $assembly -o busco -l ${params.busco_database_dir}/${params.busco_lineage} -m genome -c ${task.cpus} --offline
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

