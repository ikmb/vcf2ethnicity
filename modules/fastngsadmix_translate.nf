process FAST_NGS_ADMIX_TRANSLATE {

    publishDir "${params.outdir}/fastNGSadmix", mode: 'copy'

    tag "${meta.sample_id}"

    input:
    tuple val(meta),path(qfile)

    output:
    tuple val(meta),path(report)

    script:
    report = meta.sample_id + ".fastngsadmix.txt"

    """
        translate.rb -i $qfile > $report
    """

}
