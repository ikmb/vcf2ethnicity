process FASTNGSADMIX2EXCEL {

    container 'ikmb/exome-seq:5.4'

    publishDir "${params.outdir}/fastNGSadmix", mode: 'copy'

    tag "ALL"

    input:
    tuple val(population),path(reports)

    output:
    path(excel)

    script:
    excel = params.run_name + "_" + population + ".fastngsadmix.xlsx"

    """
        fastngsadmix2xls.rb -o $excel
    """
}