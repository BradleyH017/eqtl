# Bradley November 2023
# Filtration of genotypes for the cis-eQTL hits for the nextflow trans-by-cis module
# This will be a temporary option, which can be implemented depending on whether 

# Select vcf [This will be the input vcf for the entire workflow
repodir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl"
vcf="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/imputed_genotypes/CCF_OTAR-plates_1_2_3-imputed-all_chr.vcf.gz"

# Specify cis-eQTL results files [These are defined as the output of tensorqtl_analysis.py] - Will want to run this PER CONDITION - this does have some implications in terms of multiple testing correction
# This file has the top cis-eQTL per gene, per condition
condition="Secretory"
cis_file_path="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123" # Also will be defined within the workflow
cis_file="${cis_file_path}/category__machine-${condition}-dMean/Cis_eqtls_qval.tsv"

# Subset the variant file for those that appear in the cis-file and add to some directory
outdir=${repodir}/temp_trans
mkdir -p $outdir
# Extract variants
{ bcftools view -h $vcf; grep -Fwf <(awk '{print $7}' $cis_file) <(bcftools view $vcf);} > ${outdir}/${condition}_tmpfile.vcf

# These also need to be plink converted (Think the plink conversion module will help with this)
file__vcf=${outdir}/${condition}_tmpfile.vcf
pgen_or_bed="dosage=DS --make-pgen"
plink2_filters='--allow-extra-chr 0 --chr 1-22 XY --output-chr chrM --snps-only --rm-dup exclude-all'
plink2 --vcf ${file__vcf} ${pgen_or_bed} ${plink2_filters} --hwe 0.0000001 --out ${outdir}/${condition}_plink_genotypes

###### Also grabbed an expression bed to test this on. NOt sure if this is the correct but but doesnlt matter
${tensor_qtl_script} --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed Expression_Data.sorted.bed --covariates_file ${covariates_tsv} -window ${params.windowSize} ${dosage} --outdir ${outpath}
