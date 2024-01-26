# Bradley Jan 2024

# Load modules to run
script_path=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl/temp_scripts
export singularity=ISG/singularity/3.6.4
module load ISG/singularity/3.6.4
singularity exec -B /lustre -B /software /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl_15_05_2023.img ${script_path}/tensor_trans_by_cis_analyse.py \
    --covariates_file "Covariates.tsv" \
    --window "1000000" \
    --expression_bed "Expression_Data.bed.gz" \
    --plink_prefix_path "Secretory_plink_genotypes" \
    --outdir '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/eqtl/temp_trans' \
    --dosage \
    --maf "0.05"
