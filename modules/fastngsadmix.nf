process FAST_NGS_ADMIX {

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/fastNGSadmix", mode: 'copy'

    input:
    tuple val(meta),path(bim),path(bed),path(fam),val(ref_name),path(ref),path(ref_freq)

    output:
    tuple val(meta),val(ref_name),path(qfile), emit: results
    path("versions.yml"), emit: versions

    script:
    base_name = meta.sample_id + "-" + ref_name
    qfile = base_name + ".qopt"

    """
    fastNGSadmix -plink ${meta.sample_id} -fname $ref -Nname $ref_freq -out ${base_name} -whichPops all

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastNGSadmix: 1.0
    END_VERSIONS
    """
}
