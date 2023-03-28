process ADMIXTURE {

    label 'medium_parallel'

    publishDir "${params.outdir}/Admixture", mode: 'copy'

    tag "${meta.sample_id}"

    container 'quay.io/biocontainers/admixture:1.3.0--0'

    input:
    tuple val(meta),path(bim),path(bed),path(fam)

    output:
    tuple val(meta),path(q),path(p), emit: admix
    path("versions.yml"), emit: versions

    script:
    q = bed.getBaseName() + ".6.Q"
    p = bed.getBaseName() + ".6.P"

    """
    admixture $bed 6 -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture 1.3.0
    END_VERSIONS

    """	

}


