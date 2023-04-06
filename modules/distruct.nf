process DISTRUCT {

    label 'short_serial'

    publishDir "${params.outdir}/PDF", mode: 'copy'

    tag "${meta.sample_id}"

    input:
    tuple val(meta),path(bim),path(bed),path(fam),path(q),path(p)
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


