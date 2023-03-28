process DISTRUCT {

    label 'short_serial'

    container 'quay.io/biocontainers/pbiotools:4.0.1--pyh7cba7a3_1'

    publishDir "${params.outdir}/PDF", mode: 'copy'

    tag "${meta.sample_id}"

    input:
    tuple val(meta),path(bim),path(bed),path(fam)
    path(pops)

    output:
    tuple val(meta),path(report), emit: pdf

    script:
    title = "${meta.sample_id} - K=6"
    report = meta.sample_id + ".pdf"

    """
    distruct.py --input ${meta.sample_id} -K 6 --popfile $pops --title $title --out $report

    """	

}


