#!/usr/bin/env nextflow 
/*
    Authors:
    - Mike Sauria <msauria1@jh.edu>
*/

nextflow.enable.dsl=2
// NXF_VER=23.0" Require later version of nextflow
//assert System.getenv("NXF_VER") >= "23.0"

date = new Date().format( 'yyyyMMdd' )
contigs = Channel.from("I","II","III","IV","V","X")

if (params.debug) {
    params.sample_sheet = "${workflow.projectDir}/test_data/isotype_groups.tsv"
    params.bam_folder = "${workflow.projectDir}/test_data/bam"
    if (params.output != null) {
        params.outdir = "${params.output}"
    } else {
        params.outdir = "delly-${date}-debug"
    }
    params.species = "c_elegans"
} else if (params.help) {
    params.species = "c_elegans"
} else {
    params.species = null
}

if (params.species == "c_elegans" | params.species == "c_briggsae" | params.species == "c_tropicalis") {
    // set project and build defaults for CE, CB, and CT, can always change with argument.
    if(params.species == "c_elegans") {
        params.project="PRJNA13758"
        params.ws_build="WS283"
    } else if(params.species == "c_briggsae") {
        params.project="QX1410_nanopore"
        params.ws_build="Feb2020"
    } else if(params.species == "c_tropicalis") {
        params.project="NIC58_nanopore"
        params.ws_build="June2021"
    }
    // Define the genome
    params.genome = "${params.dataDir}/${params.species}/genomes/${params.project}/${params.ws_build}/${params.species}.${params.project}.${params.ws_build}.genome.fa.gz"
    params.genome_index = "${params.genome}.fai"
    if (params.bam_dir == null) {
        // Define bam folder
        params.bam_folder = "${params.dataDir}/data/${params.species}/WI/alignments"
    } else {
        params.bam_folder = "${params.bam_dir}"
    }
} else {
    if (params.reference == null | params.reference_index == null | params.bam_dir == null){
        params.genome = "${params.reference}"
        params.genome_index = "${params.genome}.fai"
        params.bam_folder = "${params.bam_dir}"
    } else {
        println """

            Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, or a ref genome with --reference, ref genome index with --ref_index, and bam folder with --bam_dir"

            """
            exit 1
    }
}

if (params.output != null) {
    params.outdir = "${params.output}"
} else {
    // Define output folder
    params.outdir = "delly-${date}"
}

def log_summary() {
/*
    Generates a log
*/

out = """

      ______   _______  _        _                   _        _______ 
     (  __  \\ (  ____ \\( \\      ( \\   |\\     /|     ( (    /|(  ____ \\
     | (  \\  )| (    \\/| (      | (   ( \\   / )     |  \\  ( || (    \\/
     | |   ) || (__    | |      | |    \\ (_) /_____ |   \\ | || (__    
     | |   | ||  __)   | |      | |     \\   /(_____)| (\\ \\) ||  __)   
     | |   ) || (      | |      | |      ) (        | | \\   || (      
     | (__/  )| (____/\\| (____/\\| (____/\\| |        | )  \\  || )      
     (______/ (_______/(_______/(_______/\\_/        |/    )_)|/       
                                                                                                                                                                                                                

nextflow main.nf --help

nextflow main.nf -profile rockfish --debug

nextflow main.nf -profile rockfish --sample_sheet=isotype_groups.tsv --species=c_elegans 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --sample_sheet       TSV with column isotype (needs header)                   ${params.sample_sheet}
    --output             Output folder name (optional)                            ${params.outdir}
    
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    or
    --bam_dir            Path to folder containing bams                           ${params.bam_dir}
    --reference          Path to reference genome fasta                           ${params.reference}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: https://andersenlab.org/dry-guide/latest/pipelines/pipeline-delly   
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId] 
"""
}


log.info(log_summary())


if (params.help) {
    exit 1
}


workflow { 

    // pull isotypes from sample sheet
    bams = Channel.fromPath("${params.sample_sheet}") | get_isotypes

    bams.splitText()
        .map { ["${params.bam_folder}/${it}.bam", "${params.bam_folder}/${it}.bam.bai"] }
        .combine(Channel.fromPath(params.genome))
        .combine(Channel.fromPath(params.genome_index)) | delly_call_indel
}

process get_isotypes {

    label 'xs'
    executor 'local'
    container null

    input:
        file sample_sheet

    output:
        stdout

    """
    HEADER=(`head -n 1 ${sample_sheet} | sed "s/\t/ /g"`)
    for I in \$(seq 0 1 \$((\${#HEADER[*]} - 1))); do
        if [[ \${HEADER[\${I}]} =~ ^isotype\$ ]]; then
            tail -n +2 ${sample_sheet} | cut -f \$((\${I} + 1)) | sort -k1,1 | uniq
        fi
    done
    """
}

process delly_call_indel {

    label 'md'

    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple file(bam), file(bam_index), file(genome), file(genome_index), val(out_prefix)

    output:
        file("${out_prefix}.vcf")

    """
    delly call -t INS -g ${genome} ${bam} > INS.vcf
    delly call -t DEL -g ${genome} ${bam} > DEL.vcf
    delly merge -o INDEL.vcf INS.vcf DEL.vcf
    delly filter -f somatic -o ${out_prefix}.vcf -a 0.25 INDEL.vcf
    """

}

workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    { Parameters }
    ---------------------------
    Species: ${params.species}
    bam_folder: ${params.bam_folder}
    sample_sheet: ${params.sample_sheet}
    output: ${params.outdir}
    """

    println summary

}


