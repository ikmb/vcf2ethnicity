process BCFTOOLS_MERGE {

    publishDir "${params.outdir}/Merged", mode: 'copy'

    label 'medium_parallel'

    tag "${meta.sample_id}"

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    publishDir "${params.outdir}/Bcftools", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)
    tuple path(g1k_vcf),path(vcf_g1k_tbi)

    output:
    tuple val(meta),path(merged_vcf),path(merged_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    merged_vcf = meta.sample_id + "-g1k.vcf.gz"
    merged_tbi = merged_vcf + ".tbi"
    
    """
    bcftools merge --threads ${task.cpus} -m all ${g1k_vcf} $vcf \
        | bcftools annotate -x ^FORMAT/GT,INFO \
        | bcftools view -g ^miss \
        -m2 \
        -M2 \
        -v snps \
        --min-af 0.05:minor \
        -o ${merged_vcf} \
        -O z

    bcftools index -t $merged_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}

