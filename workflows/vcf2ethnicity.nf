include { BCFTOOLS_MERGE } from '../modules/bcftools/merge'
include { PLINK_VCF ; PLINK_VCF as PLINK_VCF_COMBINED } from '../modules/plink/vcf'
include { MULTIQC } from './../modules/multiqc'
include { ADMIXTURE } from "./../modules/admixture"
include { DISTRUCT } from "./../modules/distruct"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions/main'
include { FAST_NGS_ADMIX } from "./../modules/fastngsadmix"
include { FAST_NGS_ADMIX_TRANSLATE } from "./../modules/fastngsadmix_translate"

multiqc_files = Channel.from([])
ch_versions = Channel.from([])

// ******************
// Input files and references
// *******************

ch_g1k = Channel.from( [ file("${params.g1k}", checkIfExists: true) , file("${params.g1k}.tbi", checkIfExists: true) ] )

ch_pops = Channel.from(file(params.pops,checkIfExists: true))

ch_admix_ref = Channel.from([ file(params.admix_ref), file(params.admix_freq) ])

Channel.fromPath(params.vcfs).map { v ->
	tuple( [ sample_id: file(v).getSimpleName() ], file(v, checkIfExists: true), file("${v}.tbi", checkIfExists: true ))
}.set { ch_vcfs }

// ************************************
// List of tools to run
// ************************************

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

// *******************
// The workflow starts here
// *******************

log.info "----------------------"
log.info "VCF2Ethnicity - v${workflow.manifest.version}"
log.info "----------------------"

workflow VCF2ETHNICITY {

	main:

	if ( 'fastngsadmix' in tools ) {

        PLINK_VCF(
            ch_vcfs
		)

		ch_versions = ch_versions.mix(PLINK_VCF.out.versions)
	
        FAST_NGS_ADMIX(
            PLINK_VCF.out.plink,
	    ch_admix_ref.collect()
        )

       FAST_NGS_ADMIX_TRANSLATE(
           FAST_NGS_ADMIX.out.results
       )

       ch_versions = ch_versions.mix(FAST_NGS_ADMIX.out.versions)

    } 

    if ('admixture' in tools ) {
        BCFTOOLS_MERGE(
            ch_vcfs,
            ch_g1k.collect()
        )

        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

        PLINK_VCF_COMBINED(
            BCFTOOLS_MERGE.out.vcf
        )

        ch_versions = ch_versions.mix(PLINK_VCF.out.versions)
	
        ADMIXTURE(
            PLINK_VCF.out.plink
        )

        ch_versions = ch_versions.mix(ADMIXTURE.out.versions)

        DISTRUCT(
            ADMIXTURE.out.admix,
            ch_pops.collect()
        )

	}

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    MULTIQC(
        multiqc_files.collect()
    )

	emit:
	qc = MULTIQC.out.report
	
}
