process BCFTOOLS_BCF_TO_VCF {

    label 'bcftools_bcf_to_vcf'

    input:
        tuple val(meta), path(bcf), path(bcf_index)

    output:
        tuple val(meta), path("${meta.id}_indels.vcf.gz"), path("${meta.id}_indels.vcf.gz.tbi"), emit: vcf
        path "versions.yml",                                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view -v indels -Oz5 -o ${meta.id}_indels.vcf.gz ${bcf}
    tabix -p vcf ${meta.id}_indels.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_indels.vcf.gz
    touch ${meta.id}_indels_vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}