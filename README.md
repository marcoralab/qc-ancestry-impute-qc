# QC-Ancestry-Imputation-QC *(QAIQ)* Pipeline

This pipeline is designed to be integrated as a module into custom pipelines for processing microarray gentotyping datasets. In the parent pipeline, datasets should be converted to PLINK files and probes should be filtered to select the best probes. This pipeline then does the following:

0. **Liftover**: Lift over all datasets to the same build for QC and imputation. Uses the [vcf_liftover](https://github.com/marcoralab/vcf_liftover/) pipeline.
1. **Initial sample-level QC**: Uses the [GWASampleFiltering](https://github.com/marcoralab/GWASampleFiltering) pipeline as a module in order to mark samples with low call rates, sex mismatches (based on X chromosome heterozygisity) or high overall heterozygosity. Also generates rough ancestry estimates using 1000 Genomes as a reference for PCA, then assigning continental ancestry using the geometric median of the reference populations.
2. **Ancestry Splits**: Splits the dataset into the different ancestry groups for imputation based on the results of the initial QC. This is a checkpoint and the pipeline will only submit a dataset-ancestry pair to the next step if there are samples assigned.
3. **Imputation**: Imputes the dataset using the NIH or Michigan Imputation Server. This uses [ImputePipeline]((https://github.com/marcoralab/imputePipeline)) as a module and also includes post-imputation QC.
4. **Final Sample-level QC**: Uses the [GWASampleFiltering](https://github.com/marcoralab/GWASampleFiltering) pipeline as a module to assign final PCA-based ancestry estimates, generate stratification PCs, and test for relatedness.

## Inputs from parent pipeline

Pass the config object with the `'datasets'` containing dictionary of dictionaries to define each dataset. The outer key should be the dataset name, with nested dictionaries containing the following ***keys*** and *values*:
* ***filestem***: The filestem of the PLINK files for the dataset. So if the PLINK files are named `dataset.bed`, `dataset.bim`, and `dataset.fam`, the filestem would be `dataset`.
* ***build***: The build of the dataset. Either `b37` or `b38` is prefered, but can be `hg19`, `GRCh37`, or `GRCh38` as well.
* ***nosex** (optional)*: A boolean indicating whether sex information is missing for the dataset. `False` if sex information is present, `True` if it is missing.

<br>

The code for using the QAIQ module should look like this:

```python
config_qaiq['datasets'] = {
    k: {'filestem': f'results/raw/{v}', 'build': builds[v],
        'nosex': nosex[v]}
    for k, v in zip(center_platform_dash, center_platform)}

sfile_qaic = 'modules/qc-ancestry-impute-qc/Snakefile'

module qaiq:
    snakefile: sfile_qaic
    config: config_qaiq

use rule * from qaiq as qaiq_*
```

## Other configuration options

### Requred configuration options

The pipeline requires the following additional configuration options:

* ***populations***: A list of continental ancestries to keep for imputation. Valid options include `["AFR", "AMR", "EUR", "EAS", "SAS"]`
* ***build_preimpute***: The build for initial sample QC and for upload to the imputation server. Alternatively, set the build in `Sample_Filtering_preimpute['genome_build']` for the same effect. If both are set and are different, there will be an error. This will also be used for imputation unless overridden in the `impute` configuration.
* ***chroms_preimpute***: Chromosomes for initial liftover and QC. If a dataset is `nosex`, the X chromosome will be removed for that dataset.
* ***Sample_Filtering_preimpute***: Settings for [GWASampleFiltering](https://github.com/marcoralab/GWASampleFiltering) for pre-imputation QC. `genome_build` is optional.
* ***Sample_Filtering_postimpute***: Settings for [GWASampleFiltering](https://github.com/marcoralab/GWASampleFiltering) for post-imputation QC. `genome_build` is required and should match the expected imputation server output build.
* ***impute***: The configuration options from [the ImputePipeline module](https://github.com/marcoralab/imputePipeline). In general, provide the standard configuration options, with the following differences:
  * `ref` is not required and should be left unset.
  * `chroms` will set the chromosomes for steps prior to imputation if `chroms_preimpute` is not set at the highest config level. If the X chromosome is missing form any dataset, this should not currently include `X`. This restriction will be lifted in the future.
  * `imputation`: Imputation server options. `build` is set based on build in `config['Sample_Filtering_preimpute']['genome_build']` and will be overwritten if set here. However, `token` must be provided.

  ```yaml
  impute:
    imputation:
      default: 
        token: IMPUTATION_TOKEN
        refpanel: topmed-r3
      COHORT_ANCESTRY:
        token: OTHER_IMPUTATION_TOKEN
  ```

### Optional configuration options

* ***pipeline_versions***: The paths or versions of each of the substituent modules.
* ***nointernet***: do not attempt to download references.
