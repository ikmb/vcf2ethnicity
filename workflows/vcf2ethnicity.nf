include { BCFTOOLS_MERGE } from '../modules/bcftools/merge'
include { PLINK_VCF } from '../modules/plink/vcf'
include { MULTIQC } from './../modules/multiqc'
include { ADMIXTURE } from "./../modules/admixture"
include { DISTRUCT } from "./../modules/distruct"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions/main'

multiqc_files = Channel.from([])
ch_versions = Channel.from([])

// ******************
// Input files and references
// *******************

ch_g1k = Channel.from( [ file("${params.g1k}", checkIfExists: true) , file("${params.g1k}.tbi", checkIfExists: true) ] )

ch_pops = Channel.from(file(params.pops,checkIfExists: true))

Channel.fromPath(params.vcfs).map { v ->
	tuple( [ sample_id: file(v).getSimpleName() ], file(v, checkIfExists: true), file("${v}.tbi", checkIfExists: true ))
}.set { ch_vcfs }

// *******************
// The workflow starts here
// *******************

log.info "----------------------"
log.info "VCF2Ethnicity - v${workflow.manifest.version}"
log.info "----------------------"

workflow VCF2ETHNICITY {

	main:
	
	BCFTOOLS_MERGE(
		ch_vcfs,
		ch_g1k.collect()
	)

	ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

	PLINK_VCF(
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
