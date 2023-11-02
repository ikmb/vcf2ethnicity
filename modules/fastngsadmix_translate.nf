process FAST_NGS_ADMIX_TRANSLATE {

    publishDir "${params.outdir}/fastNGSadmix/reports", mode: 'copy'

    tag "${meta.sample_id}"

    input:
    tuple val(meta),val(population),path(qfile)

    output:
    tuple val(meta),val(population),path(r), emit: report

    script:
    r = meta.sample_id + "." + population + ".fastngsadmix.txt"

    """
        translate.rb -i $qfile > $r
    """

}
