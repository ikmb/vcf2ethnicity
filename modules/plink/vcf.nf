process PLINK_VCF {

    publishDir "${params.outdir}/Plink", mode: 'copy'

    label 'medium_serial'

    tag "${meta.sample_id}"

    container 'quay.io/biocontainers/plink:1.90b6.21--hec16e2b_3'

    input:
    tuple val(meta),path(vcf),path(tbi)

    output:
    tuple val(meta),path(bim),path(bed),path(fam), emit: plink
    path("versions.yml"), emit: versions

    script:
    base = meta.sample_id
    bed = base + ".bed"
    fam = base + ".fam"
    bim = base + ".bim"

    """
    plink --vcf $vcf --make-bed --double-id --autosome --out $base

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink \$( plink --version | sed -e "s/PLINK//g" | sed -e "s/ 64-bit.*//" )
    END_VERSIONS
    """	

}


