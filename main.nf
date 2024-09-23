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
        params.refstrain="N2"
    } else if(params.species == "c_briggsae") {
        params.project="QX1410_nanopore"
        params.ws_build="Feb2020"
        params.refstrain="QC1410"
    } else if(params.species == "c_tropicalis") {
        params.project="NIC58_nanopore"
        params.ws_build="June2021"
        params.refstrain="NIC58"
    }
    // Define the genome
    params.genome = "${params.dataDir}/${params.species}/genomes/${params.project}/${params.ws_build}/${params.species}.${params.project}.${params.ws_build}.genome.fa"
        if (params.genome =~ /.*fa\.gz$/){
            params.genome_index = "${params.genome}.gzi"
        } else {
            params.genome_index = "${params.genome}.fai"
        }
    if (params.bam_dir == null) {
        // Define bam folder
        params.bam_folder = "${params.dataDir}/${params.species}/WI/alignments"
    } else {
        params.bam_folder = "${params.bam_dir}"
    }
} else {
    if (params.reference == null | params.reference_index == null | params.bam_dir == null){
        params.genome = "${params.reference}"
        if (params.genome =~ /.*fa\.gz$/){
            params.genome_index = "${params.genome}.gzi"
        } else {
            params.genome_index = "${params.genome}.fai"
        }
        params.bam_folder = "${params.bam_dir}"
    } else {
        println """

            Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, or a ref genome with --reference, ref genome index with --ref_index, and bam folder with --bam_dir"

            """
            exit 1
    }
}

if (params.masking == null) {
    if (params.species == "c_elegans") {
        params.mask_file = "${workflow.projectDir}/test_data/c_elegans_mask.bed"
    } else {
        params.mask_file = null
    }
} else {
    params.mask_file = null
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
    --masking            BED file containing regions to skip during indel calling ${params.masking}
    --minsize            The minimum size in bp to report for INDELs              ${params.minsize}
    --maxsize            The maximum size in bp to report for INDELs              ${params.maxsize}
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
    bams = Channel.fromPath("${params.sample_sheet}")
        .combine(Channel.of("${params.refstrain}")) | get_isotypes

    delly_input = bams.splitText()
        .map { it.replace("\n", "") }
        .map { ["${it}", "${params.refstrain}", file("${params.bam_folder}/${it}.bam"), file("${params.bam_folder}/${it}.bam.bai"), file("${params.bam_folder}/${params.refstrain}.bam"), file("${params.bam_folder}/${params.refstrain}.bam.bai")] }
        .combine(Channel.fromPath(params.genome))
        .combine(Channel.fromPath(params.genome_index))
        .combine(Channel.fromPath(params.mask_file))
        
    delly_call_indel( delly_input )

    delly_input
        .join(delly_call_indel.out) | delly_filter_indel
    
    delly_filter_indel.out | convert_to_vcf
}

process get_isotypes {

    label 'xs'
    executor 'local'
    container null

    input:
        tuple file(sample_sheet), val(refstrain)

    output:
        stdout

    """
    HEADER=(`head -n 1 ${sample_sheet} | sed "s/\t/ /g"`)
    for I in \$(seq 0 1 \$((\${#HEADER[*]} - 1))); do
        if [[ \${HEADER[\${I}]} =~ ^isotype\$ ]]; then
            tail -n +2 ${sample_sheet} | cut -f \$((\${I} + 1)) | sort -k1,1 | uniq | grep -v ${refstrain}
        fi
    done
    """
}

process delly_call_indel {

    label 'md'
    label "delly"

    input:
        tuple val(sample), val(control), file(bam), file(bam_index), file(ref_bam), file(ref_bam_index), file(genome), file(genome_index), file(mask)

    output:
        tuple val(sample), file("${sample}_indels_unfiltered.bcf"), file("${sample}_indels_unfiltered.bcf.csi")

    """
    if [[ -e ${mask} ]]; then MASK="-x ${mask}"; else MASK=""; fi
    delly call -q 20 -g ${genome} -o ${sample}_indels_unfiltered.bcf \$MASK ${bam} ${ref_bam}
    """
}


process delly_filter_indel {

    label 'sm'
    label "delly"

    input:
        tuple val(sample), val(control), file(bam), file(bam_index), file(ref_bam), file(ref_bam_index), file(genome), file(genome_index), file(mask), file(indels_bcf), file(indels_bcf_index)

    output:
        tuple val(sample), file("${indels_bcf}"), file("${indels_bcf_index}"), file("${sample}_indels_filtered.bcf"), file("${sample}_indels_filtered.bcf.csi")

    """
    echo -e "${control}\tcontrol\n${sample}\ttumor" > samples.tsv
    delly filter -f somatic -o ${sample}_indels_filtered.bcf -a 0.75 -p -m ${params.minsize} -n ${params.maxsize} -s samples.tsv ${indels_bcf}
    """
}

process convert_to_vcf {

    label 'xs'
    label 'annotation'

    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(sample), file(unfiltered_bcf), file(unfiltered_bcf_index), file(filtered_bcf), file(filtered_bcf_index)

    output:
        path("${sample}_*filtered_indels.vcf.gz*")

    """
    bcftools view -v indels -Oz5 -o ${sample}_unfiltered_indels.vcf.gz ${unfiltered_bcf}
    tabix -p vcf ${sample}_unfiltered_indels.vcf.gz
    bcftools view -v indels -Oz5 -o ${sample}_filtered_indels.vcf.gz ${filtered_bcf}
    tabix -p vcf ${sample}_filtered_indels.vcf.gz
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
    mask: ${params.mask_file}
    output: ${params.outdir}
    """

    println summary

}



