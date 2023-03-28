# Usage information

## Basic execution

This pipeline is configured to run on the IKMB MedCluster.

To execute, do:

```
nextflow run ikmb/vcf2ethnicity --vcfs '/path/to/*vcf.gz'
```

## Options
### `--vcfs`
This options points to one or several VCF files (compressed bgzip). Each file is assumed to be accompanied by a matching index file (.tbi).


