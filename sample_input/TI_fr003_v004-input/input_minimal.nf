params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Jun/tobi_impute/CCF_OTAR-plates_1_2_3-imputed-all_chr.vcf.gz'
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt' //assets file that has start and end positions for the genes, this one is using hg38
    phenotype_file = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/eqtl_TI_freeze_03_04_Secretory/adata.h5ad' //this should point to h5ad file in a single cell experiments.
    // This is Secretory only for testing
    // aggregation_columns='category__machine,everything__machine' //for the scrna h5ad file define which collum to use to aggregate the data based on.
    aggregation_columns='category__machine' //for the scrna h5ad file define which collum to use to aggregate the data based on.
    gt_id_column='Corrected_genotyping_ID' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='sanger_sample_id' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}