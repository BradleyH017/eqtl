tensor_label = params.TensorQTL.utilise_gpu ? 'gpu' : "process_medium"   

process TENSORQTL {  
    label "${tensor_label}"
    tag "$condition, $nr_phenotype_pcs"
    
    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}/",
                overwrite: "true"
  

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.eqtl_container}"
  } else {
    container "${params.eqtl_docker}"
  }
  

  input:
    tuple(val(condition),path(aggrnorm_counts_bed),path(covariates_tsv),val(nr_phenotype_pcs))
    each path(plink_files_prefix)

  output:
    tuple val(condition), path("${outpath}"), emit: pc_qtls_path
    path("${outpath}/Cis_eqtls_qval.tsv"), emit: qval_file // Emit the qval_file as output   
    path("${outpath}/phenotype_df.tsv"), emit: phenotype_file
    path("${outpath}/phenotype_pos_df.tsv"), emit: phenotype_pos_file


  script:
  // If a file with interaction terms is provided, use the interaction script otherwise use the standard script   
  if (params.TensorQTL.interaction_file?.trim()) {
    tensor_qtl_script = "tensorqtl_analyse_interaction.py -inter ${params.TensorQTL.interaction_file} --interaction_maf ${params.TensorQTL.interaction_maf}"
    inter_name = file(params.TensorQTL.interaction_file).baseName
    outpath = "${nr_phenotype_pcs}/interaction_output/${inter_name}"
  } else {
    tensor_qtl_script = "tensorqtl_analyse.py -nperm ${params.numberOfPermutations}"
    outpath = "${nr_phenotype_pcs}/base_output/base"
  }
  if (params.TensorQTL.use_gt_dosage) {
    dosage = "--dosage"
  }else{
    dosage = ""
  }
    """

      bedtools sort -i ${aggrnorm_counts_bed} -header > Expression_Data.sorted.bed
      ${tensor_qtl_script} --plink_prefix_path ${plink_files_prefix}/plink_genotypes --expression_bed Expression_Data.sorted.bed --covariates_file ${covariates_tsv} -window ${params.windowSize} ${dosage} --outdir ${outpath}
    """
}

// PREP_OPTIMISE_PCS process to create symlinks
process PREP_OPTIMISE_PCS {
    label 'process_low'
    tag { condition }
    input:
    tuple val(condition), val(paths)

    output:
    tuple val(condition), path("${condition}_symlink")

    script:
    paths_str = paths.join(" ")
    """
    mkdir ${condition}_symlink
    cd ${condition}_symlink
    for path in ${paths_str}; do
         unique_name=\$(echo \$path | awk -F/ '{print \$(NF-2)"__"\$(NF-1)"__"\$NF}')
        ln -s \$path \$unique_name
    done
    """
}

process OPTIMISE_PCS{
     
    // Choose the best eQTL results based on most eGenes found over a number of PCs
    // ------------------------------------------------------------------------
    tag { condition }
    scratch false      // use tmp directory
    label 'process_low'
    errorStrategy 'ignore'


    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                mode: "copy",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
        tuple(val(condition),path(eqtl_dir))
        
    output:
        tuple(
          val(condition),
          path("${outpath}/Cis_eqtls.tsv"),
          path("${outpath}/Covariates.tsv"),
          path("${outpath}/phenotype_df.tsv"),
          path("${outpath}/phenotype_pos_df.tsv"),
          path("${outpath}/Cis_eqtls_qval.tsv"),
          path("${outpath}/optimise_nPCs-FDR${alpha_text}.pdf"),
          path("${outpath}/optimise_nPCs-FDR${alpha_text}.txt"),
          path(outpath),
          emit: combined_input
        )

    script:
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
        if (params.TensorQTL.interaction_file?.trim()) {
          inter_name = file(params.TensorQTL.interaction_file).baseName
          outpath_end = "interaction_output__${inter_name}"
      } else {
          inter_name = "NA"
          outpath_end = "base_output__base"
          }
        alpha = "0.05"
        alpha_text = alpha.replaceAll("\\.", "pt")
        outpath = "./OPTIM_pcs/${outpath_end}"
        """  
          mkdir -p ${outpath}
          tensorqtl_optimise_pcs.R ./ ${alpha} ${inter_name} ${condition} ${outpath}
          var=\$(grep TRUE ${outpath}/optimise_nPCs-FDR${alpha_text}.txt | cut -f 1)
          echo \${var} >> ${outpath}/optim_pcs.txt
          cp -r ${condition}_symlink/"\$var"pcs__${outpath_end}/* ${outpath}
        """
}

process TRANS_BY_CIS {
    label "${tensor_label}"
    tag "$condition"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
    } else {
      container "${params.eqtl_docker}"
    }

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                mode: "copy",
                overwrite: "true"

    input:
        combined_input
        each path(plink_files_prefix) 

    output:
        //path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        //path("trans-by-cis_all.tsv.gz", emit: trans_res_all, optional:true)

    script:
      // Use dosage?
      if (params.TensorQTL.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }

      // Define alpha value
      alpha = "0.05"
      
      /*// Define general outdir
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
      if (params.TensorQTL.interaction_file?.trim()) {
        inter_name = file(params.TensorQTL.interaction_file).baseName
        outpath_end = "interaction_output__${inter_name}"
      } else {
        inter_name = "NA"
        outpath_end = "base_output__base"
      }
      */
      outpath = "${launchDir}/results/TensorQTL_eQTLS/${condition}/OPTIM_pcs/base_output__base/"

      """
      tensor_analyse_trans_by_cis.py \
        --covariates_file Covariates.tsv \
        --phenotype_file phenotype_df.tsv \
        --phenotype_pos_file phenotype_pos_df.tsv \
        --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
        --outdir "./" \
        --dosage ${dosage} \
        --maf "0.05" \
        --cis_qval_results Cis_eqtls_qval.tsv \
        --alpha ${alpha} \
        --window ${params.windowSize}

      """

      
}

process TRANS_OF_CIS {
    label "process_high_memory"
    tag "$condition"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.eqtl_container}"
    } else {
      container "${params.eqtl_docker}"
    }

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                mode: "copy",
                overwrite: "true"

    input:
        combined_input
        each path(plink_files_prefix) 

    output:
        //path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
        path("trans-of-cis_all.tsv", emit: trans_res, optional: true)

    script:
      // Use dosage?
      if (params.TensorQTL.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }

      // Define alpha value
      alpha = "0.05"
      
      /*// Define general outdir
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
      if (params.TensorQTL.interaction_file?.trim()) {
        inter_name = file(params.TensorQTL.interaction_file).baseName
        outpath_end = "interaction_output__${inter_name}"
      } else {
        inter_name = "NA"
        outpath_end = "base_output__base"
      }
      */
      outpath = "${launchDir}/results/TensorQTL_eQTLS/${condition}/OPTIM_pcs/base_output__base/"

      """
      tensor_analyse_trans_of_cis.py \
        --covariates_file Covariates.tsv \
        --phenotype_file phenotype_df.tsv \
        --phenotype_pos_file phenotype_pos_df.tsv \
        --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
        --outdir "./" \
        --dosage ${dosage} \
        --maf "0.05" \
        --cis_qval_results Cis_eqtls_qval.tsv \
        --alpha ${alpha} \
        --window ${params.windowSize}

      """
      //cp trans-by-cis_bonf_fdr.tsv ${outpath}
      //"""

      
}

workflow TENSORQTL_eqtls{
    take:
        condition_bed
        plink_genotype
        
    main:
  
      TENSORQTL(
          condition_bed,
          plink_genotype
      )

      if (params.TensorQTL.optimise_pcs){
          // TENSORQTL.out.pc_qtls_path.view()
          // Make sure all input files are available before running the optimisation
          TENSORQTL.out.pc_qtls_path.collect()
          // Fix the format of the output from TENSORQTL
          prep_optim_pc_channel = TENSORQTL.out.pc_qtls_path.groupTuple().map { key, values -> [key, values.flatten()] }
          // Create symlinks to the output files
          PREP_OPTIMISE_PCS(prep_optim_pc_channel)
          // Run the optimisation to get the eQTL output with the most eGenes
          OPTIMISE_PCS(PREP_OPTIMISE_PCS.out)
          channel_dsb2 = PREP_OPTIMISE_PCS.out.combine(OPTIMISE_PCS.out.combined_input, by: 0)
          if(params.TensorQTL.trans_by_cis){
            log.info 'Running trans-by-cis analysis on optimum nPCs'
            TRANS_BY_CIS(
              channel_dsb2,
              plink_genotype
            )
          }
          if(params.TensorQTL.trans_of_cis){
            log.info 'Running trans-of-cis analysis on optimum nPCs'
            TRANS_OF_CIS(
              channel_dsb2,
              plink_genotype
            )
          } 
  }
}
