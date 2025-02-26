#!/usr/bin/env nextflow 
/*
    Authors:
    - Mike Sauria <msauria1@jh.edu>
*/

include { DELLY_CALL_INDELS   } from "./modules/delly/call_indels/main"
include { DELLY_FILTER_INDELS } from "./modules/delly/filter_indels/main"
include { BCFTOOLS_BCF_TO_VCF } from "./modules/bcftools/bcf_to_vcf/main"

nextflow.enable.dsl=2
// NXF_VER=23.0" Require later version of nextflow
//assert System.getenv("NXF_VER") >= "23.0"

date = new Date().format( 'yyyyMMdd' )

if (params.debug) {
    params.sample_sheet = "${workflow.projectDir}/test_data/isotype_groups.tsv"
    bam_dir = "${workflow.projectDir}/test_data/bam"
    species = "c_elegans"
    project = "PRJNA13758"
    ws_build = "WS283"
    ref_strain = "N2"
    reference = "${params.dataDir}/${species}/genomes/${project}/${ws_build}/${species}.${project}.${ws_build}.genome.fa.gz"
    reference_index = "${params.dataDir}/${species}/genomes/${project}/${ws_build}/${species}.${project}.${ws_build}.genome.fa.gz.gzi"
} else {
    species = params.species
    if (params.bam_dir != null){
        bam_dir = params.bam_dir
    }
    if (params.ref_strain != null){
        ref_strain = params.ref_strain
    }
    if (params.reference != null){
        genome = params.reference
        if (genome =~ /.*fa\.gz$/){
            genome_index = "${pgenome}.gzi"
        } else {
            genome_index = "${genome}.fai"
        }
    }

    if (params.bam_dir == null | params.ref_strain == null | params.reference == null){
        if (species != "c_elegans" & species != "c_briggsae" & species != "c_tropicalis") {
            if (params.help == false){
                println """

                    Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, and/or a ref genome with --reference, bam directory with --bam_dir, and reference strain name with --ref_strain"

                    """
                    exit 1
            }
        } else {
            if (params.bam_dir == null){
                bam_dir = "${params.dataDir}/${species}/WI/alignments"
            }
            if (params.reference == null){
                // set project and build defaults for CE, CB, and CT, can always change with argument.
                if(species == "c_elegans") {
                    project = "PRJNA13758"
                    ws_build = "WS283"
                } else if(species == "c_briggsae") {
                    project = "QX1410_nanopore"
                    ws_build = "Feb2020"
                } else if(species == "c_tropicalis") {
                    project = "NIC58_nanopore"
                    ws_build = "June2021"
                }
                genome = "${params.dataDir}/${species}/genomes/${project}/${ws_build}/${species}.${project}.${ws_build}.genome.fa"
                genome_index = "${params.dataDir}/${species}/genomes/${project}/${ws_build}/${species}.${project}.${ws_build}.genome.fa.fai"
            }
            if (params.ref_strain == null){
                if(species == "c_elegans") {
                    ref_strain = "N2"
                } else if (species == "c_briggsae") {
                    ref_strain = "QC1410"
                } else if (species == "c_tropicalis") {
                    ref_strain = "NIC58"
                }
            }
        }
    }
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

nextflow main.nf --debug

nextflow main.nf --sample_sheet=isotype_groups.tsv --species=c_elegans 

nextflow main.nf --sample_sheet=isotype_groups.tsv --bam_dir=/path/to/bams --reference=/path/to/reference.fa --ref_strain=reference_strain 

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --sample_sheet       TSV with isotype_ref_strain column (needs header)        ${params.sample_sheet}
    --minsize            The minimum size in bp to report for INDELs              ${params.minsize}
    --maxsize            The maximum size in bp to report for INDELs              ${params.maxsize}
    
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    and / or
    --bam_dir            Path to folder containing bams                           ${bam_dir}
    --ref_strain         Name of strain to use a reference (matches genome ref)   ${ref_strain}
    --reference          Path to reference genome fasta                           ${genome}
 
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
    main:
    ch_versions = Channel.empty()

    // Pull isotypes from sample sheet, removing reference strain from isotype reference strain list
    ch_strains = Channel.splitCsv("${params.sample_sheet}", header: true, sep: "\t")
        .map{ row -> row.isotype_ref_strain }
        .unique()
        .filter{ row -> row != ref_strain}

    // Create bam channel with paths and indices
    ch_bams = ch_strains
        .map( strain -> [[id: strain], "${bam_dir}/${strain}.bam", "${bam_dir}/${strain}.bam.bai"])

    ch_ref_bam = Channel.of([[id: ref_strain], "${bam_dir}/${ref_strain}.bam", "${bam_dir}/${ref_strain}.bam.bai"])

    // Create genome channel
    ch_genome = Channel.of([genome, genome_index])

    // Call indels
    DELLY_CALL_INDELS( ch_bams,
                       ch_ref_bam,
                       ch_genome )
    ch_versions = ch_versions.mix(DELLY_CALL_INDELS.out.versions)

    // Filter indels
    DELLY_FILTER_INDELS( DELLY_CALL_INDELS.out.indels,
                         ref_strain,
                         params.minsize,
                         params.maxsize )
    ch_versions = ch_versions.mix(DELLY_FILTER_INDELS.out.versions)

    // Convert from bcf to vcf
    BCFTOOLS_BCF_TO_VCF( DELLY_FILTER_INDELS.out.filtered )
    ch_versions = ch_versions.mix(BCFTOOLS_BCF_TO_VCF.out.versions)

    // Collate and save software versions
    ch_versions
        .collectFile(name: 'workflow_software_versions.txt', sort: true, newLine: true)
        .set { ch_collated_versions }

    publish:
    BCFTOOLS_BCF_TO_VCF.out.vcf >> "indels"
    ch_collated_versions        >> "."
}

// Current bug that publish doesn't work without an output closure
output {
    "." {
        mode "copy"
    }
    "indels" {
        mode "copy"
    }
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
    Species: ${species}
    bam_dir: ${bam_dir}
    sample_sheet: ${params.sample_sheet}
    """

    println summary

}



