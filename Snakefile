# Pipeline for imputation in separate ancestries and variant QC

import re
import numpy as np
import pandas as pd
import os

os.environ["snakemake_internet"] = '1'

################################
##                            ##
##      Code for inputs       ##
##                            ##
################################

configfile: 'config/config.yaml'

version_GWASampleFiltering = config['pipeline_versions']['GWASampleFiltering']
version_imputePipeline = config['pipeline_versions']['imputePipeline']

identifiers = [x for x in config['datasets'].keys()]

def parse_chrom(chrs):
    clist = [x.split(':') for x in chrs.split(',')]
    parsed = []
    for chrs in clist:
        if len(chrs) == 2:
            chrs = [str(c) for c in range(int(chrs[0]), int(chrs[1]) + 1)]
        elif len(chrs) != 1:
            raise ValueError('Invalid chromosome list.')
        parsed += chrs
    return parsed


chromlist_preimpute = parse_chrom(config['impute']['chroms'])

def flatten(nested):
    flat = []
    for el in nested:
        if not isinstance(el, list):
            flat.append(el)
        else:
            flat += flatten(el)
    return flat

imputed_stem = 'intermediate/imputation/imputed/processed'
imputed_outs = dict(
    stat_report='results/imputed/stats/{cohort}_imputation.html',
    vcf_bycohort='results/imputed/{cohort}.vcf.gz',
    vcf_merged='results/imputed/all.vcf.gz',
    bgen_bycohort='{imputed_stem}/data/{cohort}_chrall_filtered.bgen',
    bgen_merged='{imputed_stem}/data/merged/merged_chrall_filtered.bgen',
    plink_bycohort='results/imputed/{cohort}.{ext}',
    plink_merged='results/imputed/all.{ext}')

identifier_ancestry = expand('{identifier}_{ancestry}',
                             imputed_stem=imputed_stem,
                             identifier=identifiers,
                             ancestry=config['populations'])

def identifier_ancestry_present(wc):
    cp = checkpoints.count_spop.get(**wc)
    with cp.output[0].open() as f:
        count_tab = pd.read_table(f)
    keep_row = [x in config['populations']
                for x in count_tab['ancestry']]
    count_tab_filt = count_tab[keep_row]
    return count_tab_filt['identifier_ancestry'].tolist()

def imputed_outputs(wc):
    def expand_outs(out, id_a):
        return expand(out, cohort=id_a,
                      ext=['bed', 'bim', 'fam'])
    id_a = identifier_ancestry_present(wc)
    return flatten([expand_outs(imputed_outs[x], id_a)
                    for x in config["impute"]["outputs"]])

################################
##                            ##
##      Snakemake setup       ##
##                            ##
################################

shell.executable('/bin/bash')

wildcard_constraints:
    identifier = r'[^/]+'

localrules: all, imputation_submit_imputation, imputation_download_imputation, ref_download_md5_b38, ref_download_md5_hg19, ref_download_tg_fa, ref_download_tg_ped, ref_download_tg_chrom, ref_download_md5_b38, ref_download_md5_hg19

rule all:
    input:
        expand(
            'results/preimpute_stats/{identifier}_qc.html',
            identifier=identifiers),
        imputed_outputs,
        'results/imputed/stats/all_qc.html',
        'results/exclude.samples'
    output:
        touch('results/qaic.done')


########################################
###                                  ###
###     Pre-imputation Sample QC     ###
###                                  ###
########################################

## GWAS Sample Filtering References
refs = {}
refs['ref_only'] = True
refs['genome_build'] = ['hg19', 'hg38']
refs['nointernet'] = True #TODO Delete this
refs['GenoMiss'] = list(set([
    config['Sample_Filtering_preimpute']['QC']['GenoMiss'],
    config['Sample_Filtering_postimpute']['QC']['GenoMiss']]))

local_src = '/sc/arion/projects/load/users/fultob01/src'

if version_GWASampleFiltering == 'local':
    sfile_GWASampleFiltering_ref = local_src + '/GWASampleFiltering/workflow/rules/reference.smk'
elif re.search(r'/', version_GWASampleFiltering):
    sfile_GWASampleFiltering_ref = os.path.join(
      os.path.abspath(version_GWASampleFiltering), 'workflow', 'rules', 'reference.smk')
else:
    sfile_GWASampleFiltering_ref = github('marcoralab/GWASampleFiltering', path='workflow/rules/reference.smk', tag=version_GWASampleFiltering)

module refs:
    snakefile: sfile_GWASampleFiltering_ref
    config: refs

use rule * from refs as ref_*

## GWAS Sample Filtering
prefilter = config['Sample_Filtering_preimpute']
prefilter['is_module'] = True
prefilter['construct_ref'] = False
prefilter['nointernet'] = True
prefilter['SAMPLE'] = identifiers
prefilter['start'] = {'files': multiext('input/{sample}',
                                        '.bed', '.bim', '.fam'),
                       'stem': 'input/{sample}'}
prefilter['start']['sex'] = prefilter['start']['files']
prefilter['start']['sex_stem'] = prefilter['start']['stem']
prefilter['DATAOUT'] = 'output'
prefilter_prefix = 'intermediate/pre-impute_filter'

if version_GWASampleFiltering == 'local':
    sfile_GWASampleFiltering = local_src + '/GWASampleFiltering/workflow/Snakefile'
elif re.search(r'/', version_GWASampleFiltering):
    sfile_GWASampleFiltering = os.path.join(
      os.path.abspath(version_GWASampleFiltering), 'workflow', 'Snakefile')
else:
    sfile_GWASampleFiltering = github('marcoralab/GWASampleFiltering', path='workflow/Snakefile', tag=version_GWASampleFiltering)

module prefiltering:
    snakefile: sfile_GWASampleFiltering
    config: prefilter
    prefix: prefilter_prefix

use rule * from prefiltering as prefiltering_*

rule link_unimputed:
    input:
        bed = lambda wc: config['datasets'][wc.identifier] + '.bed',
        bim = lambda wc: config['datasets'][wc.identifier] + '.bim',
        fam = lambda wc: config['datasets'][wc.identifier] + '.fam'
    output:
        bed = prefilter_prefix + '/input/{identifier}.bed',
        bim = prefilter_prefix + '/input/{identifier}.bim',
        fam = prefilter_prefix + '/input/{identifier}.fam'
    threads: 1
    resources:
        mem_mb = 512,
        time_min = 5
    shell:
        '''
ln -sr {input.bed} {output.bed}
ln -sr {input.bim} {output.bim}
ln -sr {input.fam} {output.fam}
'''

rule link_prefilt_refname:
    input:
        pops = 'reference/{refname}_pops.txt',
        pops_unique = 'reference/{refname}_pops_unique.txt',
    output:
        pops = prefilter_prefix + '/' + 'reference/{refname}_pops.txt',
        pops_unique = prefilter_prefix + '/' + 'reference/{refname}_pops_unique.txt',
    resources:
        time_min = 2
    shell:
        '''
ln -rs {input.pops} {output.pops}
ln -rs {input.pops_unique} {output.pops_unique}
'''

rule link_prefilt_refname_gbuild_miss:
    input:
        vcf = 'reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz',
        snps = 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
    output:
        vcf = prefilter_prefix + '/' + 'reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz',
        snps = prefilter_prefix + '/' + 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
    resources:
        time_min = 2
    shell:
        '''
ln -rs {input.vcf} {output.vcf}
ln -rs {input.snps} {output.snps}
'''

rule link_prefilt_fasta:
    input: 'reference/human_g1k_{gbuild}.fasta',
    output: prefilter_prefix + '/' + 'reference/human_g1k_{gbuild}.fasta'
    resources:
        time_min = 2
    shell: 'ln -rs {input} {output}'

rule link_prefilt_ped:
    input: 'reference/20130606_g1k.ped'
    output: prefilter_prefix + '/' + 'reference/20130606_g1k.ped'
    resources:
        time_min = 2
    shell: 'ln -rs {input} {output}'


rule ancestry_sample_lists:
    input: prefilter_prefix + '/output/{identifier}_cluster_pops.tsv'
    output: 'intermediate/pre-impute_ancestry/{identifier}_{ancestry}.samples'
    conda: 'envs/miller.yaml'
    threads: 1
    resources:
        mem_mb = 1024,
        time_min = 180
    shell:
        '''
mlr --itsv --onidx --fs tab \
  filter '$superpop_infered == "{wildcards.ancestry}"' \
  then cut -f FID,IID {input} > {output}
'''

checkpoint count_spop:
    input:
        expand('intermediate/pre-impute_filter/output/{identifier}_cluster_pops.tsv',
               identifier=identifiers)
    output: 'intermediate/pre-impute_filter/output/ancestry_counts.tsv'
    conda: 'envs/R.yaml'
    threads: 1
    resources:
        mem_mb = 4096,
        time_min = 10
    script: 'rule_count_spop.R'

rule ancestry_sample_filt:
    input:
        samps = rules.ancestry_sample_lists.output,
        plink = lambda wc: [config['datasets'][wc.identifier] + x for x in ['.bed', '.bim', '.fam']]
    output:
        multiext('intermediate/pre-impute_ancestry/{identifier}_{ancestry}',
                 '.bed', '.bim', '.fam')
    params:
        ins = lambda wc: config['datasets'][wc.identifier],
        out = 'intermediate/pre-impute_ancestry/{identifier}_{ancestry}'
    conda: 'envs/PLINK.yaml'
    threads: 1
    resources:
        mem_mb = 4096,
        time_min = 180
    shell:
        '''
plink --bfile {params.ins} --keep-allele-order --keep {input.samps} \
  --make-bed --out {params.out}
'''

########################################
###                                  ###
###   Imputation                     ###
###                                  ###
########################################

config['impute']['out_dir'] = 'intermediate/imputation/'
config['impute']['directory'] = 'intermediate/pre-impute_ancestry'
config['impute']['SAMPLES'] = identifier_ancestry # [x.replace('_', '.') for x in center_platform_ancestry]
genome_build_preimp = config['Sample_Filtering_preimpute']['genome_build'].lower()
if genome_build_preimp in ['hg19', 'hg37', 'grch37', 'b37']:
    config['impute']['ref'] = 'reference/human_g1k_hg19.fasta'
elif genome_build_preimp in ['hg38', 'grch38', 'b38']:
    config['impute']['ref'] = 'reference/human_g1k_GRCh38.fasta'
else:
    raise ValueError('genome build must be hg19 or hg38')

if 'postImpute' in config['pipeline_versions']:
    config['impute']['version_postImpute'] = config['pipeline_versions']['postImpute']
elif version_imputePipeline == 'local':
    config['impute']['version_postImpute'] = local_src + '/imputePipeline/workflow/modules/postImpute/'

if version_imputePipeline == 'local':
    sfile_imputation = local_src + '/imputePipeline/workflow/Snakefile'
elif re.search(r'/', version_imputePipeline):
    sfile_imputation = os.path.join(
      os.path.abspath(version_imputePipeline), 'workflow', 'Snakefile')
else:
    sfile_imputation = github('marcoralab/imputePipeline', path='workflow/Snakefile', tag=version_imputePipeline)


module imputation:
    snakefile: sfile_imputation
    config: config['impute']

use rule * from imputation as imputation_*

rule cat_fams:
    input: [config['datasets'][x] + '.fam' for x in identifiers]
    output: 'intermediate/all_preimpute.fam'
    threads: 1
    resources:
        mem_mb = 1024,
        time_min = 180
    shell: r'''
awk 'BEGIN {{
    OFS = "\t"
}}
NF != 6 {{
    print "Error: Number of fields on line " FNR " of " FILENAME " is not 6"
    exit 1
}}
1 {{
    $1 = $1
    print
}}' {input} > {output}
'''

rule fix_fam:
    input:
        oldfam = rules.cat_fams.output,
        newfam = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered.fam'
    output: 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered_fixed.fam'
    threads: 1
    resources:
        mem_mb = 1024,
        time_min = 180
    conda: "envs/R.yaml"
    script: 'fix_fam.R'

use rule link_unimputed as link_imputed_output with:
    input:
        bed = 'intermediate/imputation/imputed/processed/data/all_chrall_filtered.bed',
        bim = 'intermediate/imputation/imputed/processed/data/all_chrall_filtered.bim',
        fam = expand(rules.fix_fam.output, cohort="all")
    output:
        bed = 'results/imputed/all.bed',
        bim = 'results/imputed/all.bim',
        fam = 'results/imputed/all.fam'

use rule link_unimputed as link_imputed_output_cohort with:
    input:
        bed = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered.bed',
        bim = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered.bim',
        fam = rules.fix_fam.output
    output:
        bed = 'results/imputed/{cohort}.bed',
        bim = 'results/imputed/{cohort}.bim',
        fam = 'results/imputed/{cohort}.fam'

rule link_imputed_vcf:
    input: 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered.vcf.gz'
    output: 'results/imputed/{cohort}.vcf.gz'
    threads: 1
    resources:
        mem_mb = 512,
        time_min = 2
    shell:
        '''
ln -sr {input} {output}
'''

use rule link_imputed_vcf as link_imputed_stats with:
    input: 'intermediate/imputation/imputed/processed/stats/{cohort}_impStats.html'
    output: 'results/imputed/stats/{cohort}_imputation.html'

use rule link_unimputed as link_imputed with:
    input:
        bed = rules.link_imputed_output.output.bed,
        bim = rules.link_imputed_output.output.bim,
        fam = rules.link_imputed_output.output.fam
    output:
        bed = 'intermediate/post-impute_filter/input/imputed.bed',
        bim = 'intermediate/post-impute_filter/input/imputed.bim',
        fam = 'intermediate/post-impute_filter/input/imputed.fam'

use rule stats from postImpute as postImpute_stats with:
    resources:
        mem_mb = 24000,
        walltime = '8:00'

def imputation_getmerge(wc):
    # defaults for renaming:
    renamed_merge = "{{impute_dir}}/data/by_chrom/{cohort}_chr{{chrom}}_filtered.vcf.gz"
    #override defaults:
    if 'rename' in config['impute']:
        if (config['impute']['rename'] is not None and
            'automap' in config['impute']['rename'] and
            config['impute']['rename']['automap']):
            renamed_merge = "{{impute_dir}}/data/by_chrom/{cohort}_chr{{chrom}}_filtered_fixedIDs.vcf.gz"
        elif config['impute']['rename'] and not type(config['impute']['rename']) is dict:
            renamed_merge = "{{impute_dir}}/data/by_chrom/{cohort}_chr{{chrom}}_filtered_renamed.vcf.gz"
    COHORT = identifier_ancestry_present(wc)
    return expand(renamed_merge, cohort=COHORT)

def imputation_getmerge_bgen(wc):
    COHORT = identifier_ancestry_present(wc)
    renamed_merge = "{{impute_dir}}/temp/{cohort}_chr{{chrom}}"
    return expand(renamed_merge, cohort=COHORT)

use rule merge_samples_chrom from postImpute as postImpute_merge_samples_chrom with:
    input:
        vcf = lambda wc: imputation_getmerge(wc),
        tbi = lambda wc: [x + ".tbi" for x in imputation_getmerge(wc)]

if 'bgen_merged' in config['impute']['outputs']:
    use rule make_bgen_allsamp from postImpute as postImpute_make_bgen_allsamp with:
        input:
            gen = lambda wc: [x + "_filtered.bgen" for x in imputation_getmerge_bgen(wc)],
            samp = lambda wc: [x + ".sample" for x in imputation_getmerge_bgen(wc)]

########################################
###                                  ###
###     Pre-imputation Sample QC     ###
###                                  ###
########################################

## GWAS Sample Filtering
postfilter = config['Sample_Filtering_postimpute']
postfilter['is_module'] = True
postfilter['construct_ref'] = False
postfilter['SAMPLE'] = ['imputed']
postfilter['start'] = {'files': multiext('input/{sample}',
                                        '.bed', '.bim', '.fam'),
                       'stem': 'input/{sample}'}
postfilter['start']['sex'] = postfilter['start']['files']
postfilter['start']['sex_stem'] = postfilter['start']['stem']
postfilter['DATAOUT'] = 'output'
postfilter['nointernet'] = True
postfilter_prefix = 'intermediate/post-impute_filter'


module postfiltering:
    snakefile: sfile_GWASampleFiltering
    config: postfilter
    prefix: postfilter_prefix

use rule * from postfiltering as postfiltering_*


use rule Sample_Flip from postfiltering as postfiltering_Sample_Flip with:
    threads: 10
    resources:
        mem_mb = 6000,
        walltime = '8:00'

rule link_postfilt_refname:
    input:
        pops = 'reference/{refname}_pops.txt',
        pops_unique = 'reference/{refname}_pops_unique.txt',
    output:
        pops = postfilter_prefix + '/' + 'reference/{refname}_pops.txt',
        pops_unique = postfilter_prefix + '/' + 'reference/{refname}_pops_unique.txt',
    resources:
        time_min = 2
    shell:
        '''
ln -rs {input.pops} {output.pops}
ln -rs {input.pops_unique} {output.pops_unique}
'''

rule link_postfilt_refname_gbuild_miss:
    input:
        vcf = 'reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz',
        snps = 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
    output:
        vcf = postfilter_prefix + '/' + 'reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz',
        snps = postfilter_prefix + '/' + 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
    resources:
        time_min = 2
    shell:
        '''
ln -rs {input.vcf} {output.vcf}
ln -rs {input.snps} {output.snps}
'''

rule link_postfilt_fasta:
    input: 'reference/human_g1k_{gbuild}.fasta',
    output: postfilter_prefix + '/' + 'reference/human_g1k_{gbuild}.fasta'
    resources:
        time_min = 2
    shell: 'ln -rs {input} {output}'

rule link_postfilt_ped:
    input: 'reference/20130606_g1k.ped'
    output: postfilter_prefix + '/' + 'reference/20130606_g1k.ped'
    resources:
        time_min = 2
    shell: 'ln -rs {input} {output}'

use rule link_imputed_vcf as link_imputed_qcstats with:
    input: 'intermediate/post-impute_filter/output/stats/imputed_GWAS_QC.html'
    output: 'results/imputed/stats/all_qc.html'

use rule link_imputed_vcf as link_unimputed_qcstats with:
    input: 'intermediate/pre-impute_filter/output/stats/{identifier}_GWAS_QC.html'
    output: 'results/preimpute_stats/{identifier}_qc.html'

rule cat_exclusions:
    input:
        "intermediate/post-impute_filter/output/imputed_exclude.samples",
        expand(
            'intermediate/pre-impute_filter/output/{identifier}_exclude.samples',
            identifier=identifiers)
    output: 'results/exclude.samples'
    shell:
        '''
awk 'NR==1 || FNR>1' {input} > {output}
'''
