process FAST_NGS_ADMIX {

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/fastNGSadmix", mode: 'copy'

    input:
    tuple val(meta),path(bim),path(bed),path(fam)
    tuple path(ref),path(ref_freq)

    output:
    tuple val(meta),path(qfile), emit: results
    path("versions.yml"), emit: versions

    script:

    qfile = meta.sample_id + ".qopt"

    """
    /work_beegfs/ikmb_repository/software/fastngsadmix/1.0/source/fastNGSadmix -plink "${meta.sample_id} -fname $ref -Nname $ref_freq -out ${meta.sample_id} -whichPops all"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastNGSadmix 1.0
    END_VERSIONS
    """
}