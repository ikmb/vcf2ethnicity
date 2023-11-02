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

### tools [default = 'fastngsadmix']
This option allows you to specify the desired analysis tool chain(s). Valid choices are:

* [fastNGSadmix](http://www.popgen.dk/software/index.php/FastNGSadmix) (fastngsadmix)
* [Admixture](https://dalexander.github.io/admixture/) (admixture)

These can be provided like so:

```
nextflow run ikmb/vcf2ethnicity --tools 'fastngsadmix,admixture' --vcfs '/path/to/*.vcf.gz'
```

We recommend to use fastNGSadmix only, since the overall runtime is much shorter and results are likely to be equivalent between the tools. 