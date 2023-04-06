process PLINK_VCF {

    publishDir "${params.outdir}/Plink", mode: 'copy'

    label 'medium_serial'

    tag "${meta.sample_id}"

    container 'quay.io/biocontainers/plink2:2.00a3.7--h9f5acd7_2'

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
    plink2 --vcf $vcf --snps-only --make-bed --out $base --autosome --double-id --max-alleles 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$( plink2 --version | sed -e "s/PLINK//g" | sed -e "s/ 64-bit.*//" )
    END_VERSIONS
    """	

}


