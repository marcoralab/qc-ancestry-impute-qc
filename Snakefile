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
version_vcf_liftover = config['pipeline_versions']['vcf_liftover']

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

localrules: all #, imputation_submit_imputation, imputation_download_imputation, ref_download_md5_b38, ref_download_md5_hg19, ref_download_tg_fa, ref_download_tg_ped, ref_download_tg_chrom, ref_download_md5_b38, ref_download_md5_hg19

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
###     Pre-imputation Liftover      ###
###                                  ###
########################################

## Preparation steps

def liftover_align_fasta(wc):
    build = config['datasets'][wc.identifier]['build']
    if build.lower() in ['b37', 'hg19', 'grch37']:
        return f"resources/ref/b37.fa.gz"
    elif build.lower() in ['b38', 'hg38', 'grch38']:
        return f"resources/ref/b38.fa.gz"
    else:
        raise ValueError(f"Invalid build {build}")

#add Sample flip
rule liftover_align:
    input:
        bed = lambda wc: config['datasets'][wc.identifier]['filestem'] + '.bed',
        bim = lambda wc: config['datasets'][wc.identifier]['filestem'] + '.bim',
        fam = lambda wc: config['datasets'][wc.identifier]['filestem'] + '.fam',
        fasta = liftover_align_fasta
    output:
        multiext("intermediate/liftover_align/{identifier}_flipped", ".bim", ".bed", ".fam")
    params:
        outstem = "intermediate/liftover_align/{identifier}"
    resources:
        mem_mb = 10000,
        time_min = 30
    container: 'docker://befh/flippyr:0.5.3'
    shell: "flippyr -p {input.fasta} -o {params.outstem} {input.bim}"

# Recode sample plink file to vcf
rule liftover_Plink2Vcf:
    input:
        bed = "intermediate/liftover_align/{identifier}_flipped.bed",
        bim = "intermediate/liftover_align/{identifier}_flipped.bim",
        fam = "intermediate/liftover_align/{identifier}_flipped.fam"
    output: "intermediate/liftover_align/{identifier}_flipped.vcf.gz"
    params:
        outstem = "intermediate/liftover_align/{identifier}_flipped"
    resources:
        mem_mb = 10000,
        time_min = 180
    conda: "envs/PLINK.yaml"
    shell:
        '''
        plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode vcf bgz --real-ref-alleles --out {params.outstem}
	'''

def rename_chrs(wc):
    build = config['datasets'][wc.identifier]['build']
    if build.lower() in ['b37', 'hg19', 'grch37']:
        return f"reference/rename_chrs_b37.txt"
    elif build.lower() in ['b38', 'hg38', 'grch38']:
        return f"reference/rename_chrs_b38.txt"
    else:
        raise ValueError(f"Invalid build {build}")

rule rename_chrs:
    input: "intermediate/liftover_align/{identifier}_flipped.vcf.gz"
    output: "intermediate/liftover_align/{identifier}_flipped_renamed.vcf.gz"
    params:
        txt = rename_chrs
    resources:
        mem_mb = 10000,
        time_min = 120
    #conda: "envs/bcftools.yaml"
    shell:
        '''
        ml bcftools
        bcftools annotate --rename-chrs {params.txt} -Oz -o {output} {input}
        bcftools index -t {output}
	'''

## VCF liftover module
liftover = {}
liftover['inputs'] = {
    k: {'input_vcf': f"intermediate/liftover_align/{k}_flipped_renamed.vcf.gz",
        'build': config['datasets'][k]['build'],
        'contigs': config['impute']['chroms']} 
    for k in identifiers
}
liftover['output_builds'] = ['b38']
liftover['all_GATK_builds'] = True
liftover['concatenate'] = True
liftover['liftover_outdir'] = "intermediate/lifted"

if version_vcf_liftover == 'local':
    sfile_vcf_liftover = local_src + '/vcf_liftover/workflow/Snakefile'
elif re.search(r'/', version_vcf_liftover):
    sfile_vcf_liftover = os.path.join(
      os.path.abspath(version_vcf_liftover), 'workflow', 'Snakefile')
else:
    sfile_vcf_liftover = github('marcoralab/vcf_liftover', path='workflow/Snakefile', tag=version_vcf_liftover)

module vcf_liftover:
    snakefile: sfile_vcf_liftover
    config: liftover

use rule * from vcf_liftover as vcf_liftover_*

rule make_lifted_plink_all:
    input: 'intermediate/lifted/{identifier}.b38.same_chr.vcf.gz'
    output: multiext("intermediate/lifted/{identifier}", ".bed", ".bim", ".fam")
    params:
        out_plink = "intermediate/lifted/{identifier}"
    threads: 10
    resources:
        mem_mb = 3000,
        walltime = "96:00"
    conda: "envs/PLINK.yaml"
    shell: 'plink --keep-allele-order --vcf {input} --double-id --memory 10000 --threads 10 --make-bed --out {params.out_plink}'

rule cat_fams:
    input: [ config['datasets'][x]['filestem'] + '.fam' for x in identifiers ]
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
    output: 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered_fixed.famx'
    threads: 1
    resources:
        mem_mb = 1024,
        time_min = 180
    conda: "envs/R.yaml"
    script: 'fix_fam.R'

use rule fix_fam as fix_lifted_fam with:
    input: 
        oldfam = rules.cat_fams.output,
        newfam = 'intermediate/lifted/{identifier}.fam'
    output:
        "intermediate/lifted/{identifier}_fixed.fam"


########################################
###                                  ###
###     Pre-imputation Sample QC     ###
###                                  ###
########################################

prefilter_prefix = 'intermediate/pre-impute_filter'

rule link_unimputed:
    input:
        bed = "intermediate/lifted/{identifier}.bed",
        bim = "intermediate/lifted/{identifier}.bim",
        fam = "intermediate/lifted/{identifier}_fixed.fam"
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
  filter '$superpop_infered == "{wildcards.ancestry}" && $cohort != "Reference"' \
  then put '$FID_IID = format("{{}}_{{}}", $FID, $IID)' \
  then cut -f FID_IID {input} > {output}
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


########################################
###                                  ###
###   Imputation                     ###
###                                  ###
########################################

## Preparation for imputation

rule ancestry_sample_filt:
    input:
        samps = rules.ancestry_sample_lists.output,
        vcf = 'intermediate/lifted/{identifier}.b{build}.chr{chrom}_only.vcf.gz'
    output:
        vcf = temp('intermediate/pre-impute_ancestry/{identifier}_{ancestry}.b{build}.chr{chrom}.vcf.gz'),
        tbi = 'intermediate/pre-impute_ancestry/{identifier}_{ancestry}.b{build}.chr{chrom}.vcf.gz.tbi'
    conda: 'envs/bcftools.yaml'
    resources:
        mem_mb = 10000,
        time_min = 120
    shell: '''
bcftools view -S {input.samps} -Oz -o {output.vcf} {input.vcf}
bcftools index -tf {output.vcf}
'''

rule rename_for_dup:
    input: rules.ancestry_sample_filt.output.vcf
    output: temp('intermediate/formerge/{identifier}_{ancestry}_b{build}_replicate{dup}.chr{chrom}.vcf.gz')
    conda: 'envs/bcftools.yaml'
    resources:
        mem_mb = 10000,
        time_min = 120
    shell: '''
ln -rs {input} {output}
bcftools index -tf {output}
'''

rule dup_samples:
    input:
        vcf = rules.ancestry_sample_filt.output.vcf,
        vcfvcfvcf = lambda wc: expand('intermediate/formerge/{{identifier}}_{{ancestry}}_b{{build}}_replicate{dup}.chr{{chrom}}.vcf.gz',
                                      dup=[x + 1 for x in range(int(wc['rep']) - 1)]),
        tbi = rules.ancestry_sample_filt.output.tbi
    output: temp('intermediate/pre-impute_ancestry/{identifier}_{ancestry}_b{build}_duplicated{rep}.chr{chrom}.vcf.gz')
    conda: 'envs/bcftools.yaml'
    resources:
        mem_mb = 10000,
        time_min = 120
    shell: 'bcftools merge --force-samples {input.vcf} {input.vcfvcfvcf} -Oz -o {output}'


def identifier_ancestry_dup(wc):
    cp = checkpoints.count_spop.get(**wc)
    with cp.output[0].open() as f:
        count_tab = pd.read_table(f)
    count_tab['dups'] = np.int32(np.ceil(30 / count_tab['n']))
    ia = count_tab[count_tab['identifier_ancestry'] == f"{wc['identifier']}_{wc['ancestry']}"]
    assert ia.shape[0] == 1
    dups = ia.loc[ia.index[0], 'dups']
    if dups == 1:
        return expand('intermediate/pre-impute_ancestry/{identifier}_{ancestry}.b{build}.chr{chrom}.vcf.gz', **wc)
    else:
        return expand('intermediate/pre-impute_ancestry/{identifier}_{ancestry}_b{build}_duplicated{rep}.chr{chrom}.vcf.gz',
                      identifier=wc['identifier'], ancestry=wc['ancestry'], build=wc['build'], chrom=wc['chrom'], rep=dups)

rule copy_imputation_input:
    input: identifier_ancestry_dup
    output: temp('intermediate/imputation/input/{identifier}_{ancestry}.b{build}.chr{chrom}.vcf.gz')
    localrule: True
    shell: 'cp {input} {output}'


rule make_keep_lists:
    input: 'intermediate/liftover_align/{identifier}_flipped_renamed.vcf.gz' #or lifted 'intermediate/lifted/{identifier}.b38.same_chr.vcf.gz'
    output: 'intermediate/sample_lists/{identifier}_keep_samples.txt'
    localrule: True
    conda: 'envs/bcftools.yaml'
    shell: 'bcftools query --list-samples {input} > {output}'


rule cat_keep_lists:
    input: expand('intermediate/sample_lists/{identifier}_keep_samples.txt', identifier=identifiers)
    output: 'intermediate/sample_lists/keep_samples_all.txt'
    localrule: True
    shell: 'cat {input} > {output}'


## Imputation Module

config['impute']['out_dir'] = 'intermediate/imputation/'
config['impute']['directory'] = 'intermediate/imputation/input'
config['impute']['include_samp_post'] = str(rules.cat_keep_lists.output[0])

genome_build_preimp = config['Sample_Filtering_preimpute']['genome_build'].lower()

if genome_build_preimp in ['hg19', 'hg37', 'grch37', 'b37']:
    config['impute']['ref'] = 'resources/ref/b37.fa.gz'
    config['impute']['imputation']['default']['build'] = 'hg19'
    ibuild = "37"
elif genome_build_preimp in ['hg38', 'grch38', 'b38']:
    config['impute']['ref'] = 'resources/ref/b38.fa.gz'
    config['impute']['imputation']['default']['build'] = 'hg38'
    ibuild = "38"
else:
    raise ValueError('genome build must be hg19 or hg38')

config['impute']['SAMPLES'] = {
    ia: {'file': f'intermediate/imputation/input/{ia}.b{ibuild}.chr{{chrom}}.vcf.gz',
         'type': 'vcf_chr'} for ia in identifier_ancestry}

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

#rule cat_fams happens here (moved up in pipeline for use among liftover rules)

#rule fix_fam happens here (moved up in pipeline for use among liftover rules)

use rule link_unimputed as link_imputed_output with:
    input:
        bed = 'intermediate/imputation/imputed/processed/data/all_chrall_filtered_fixed.bed',
        bim = 'intermediate/imputation/imputed/processed/data/all_chrall_filtered_fixed.bim',
        fam = 'intermediate/imputation/imputed/processed/data/all_chrall_filtered_fixed.fam',
    output:
        bed = 'results/imputed/all.bed',
        bim = 'results/imputed/all.bim',
        fam = 'results/imputed/all.fam'

use rule link_unimputed as link_imputed_output_cohort with:
    input:
        bed = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered_fixed.bed',
        bim = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered_fixed.bim',
        fam = 'intermediate/imputation/imputed/processed/data/{cohort}_chrall_filtered_fixed.fam',
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
###     Post-imputation Sample QC    ###
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
