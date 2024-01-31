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
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.pdf"), emit: optimise_nPCs_plot
        path("${outpath}/optimise_nPCs-FDR${alpha_text}.txt"), emit: optimise_nPCs
        path("${outpath}/Cis_eqtls.tsv"), emit: optim_qtl_bin, optional: true
        path("${outpath}/Cis_eqtls_qval.tsv"), emit: optim_q_qtl_bin, optional: true
        path("${outpath}/cis_inter1.cis_qtl_top_assoc.txt.gz "), emit: optim_int_qtl_bin, optional: true
        path("${outpath}/optim_pcs.txt"), emit: optim_var, optional:true
        path("${outpath}/Covariates.tsv"), emit: optim_covariates, optional:true
        path(outpath), emit: out_path, optional: true
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
                mode: "copy"

    input:
        tuple(val(condition),path(aggrnorm_counts_bed))
        each path(plink_files_prefix) 
        path(optim_q_qtl_bin)
        path(optim_covariates)
        path(out_path)
    
    output:
        path("${condition}_test_params.txt")
        path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
    
    script:
      // Use dosage?
      if (params.TensorQTL.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }

      // Define alpha value
      alpha = "0.05"

      // Define general outdir
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
      if (params.TensorQTL.interaction_file?.trim()) {
        inter_name = file(params.TensorQTL.interaction_file).baseName
        outpath_end = "interaction_output__${inter_name}"
      } else {
        inter_name = "NA"
        outpath_end = "base_output__base"
      }
      outpath = "../../../${sumstats_path}OPTIM_pcs/${outpath_end}"

      """
      tensor_analyse_trans_by_cis.py \
        --covariates_file ${optim_covariates} \
        --expression_bed Expression_Data.bed.gz \
        --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
        --outdir ${out_path} \
        --dosage ${dosage} \
        --maf "0.05" \
        --cis_qval_results ${optim_q_qtl_bin} \
        --alpha ${alpha} \
        --window ${params.windowSize}

      cp -r ${out_path}/* ${outpath}
      """

      
}

/*
// Execute trans by cis
process TRANS_BY_CIS_OLD{

    // Subset for significant cis-eQTL effects and perform trans-by-cis analysis within each condition
    // ------------------------------------------------------------------------
    label "${tensor_label}"
    tag "$condition, $nr_phenotype_pcs"

    publishDir  path: "${params.outdir}/TensorQTL_eQTLS/${condition}",
                mode: "copy",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.eqtl_container}"
    } else {
        container "${params.eqtl_docker}"
    }

    input:
      tuple(val(condition),path(aggrnorm_counts_bed),path(covariates_tsv),val(nr_phenotype_pcs))
      each path(plink_files_prefix)
      tuple val(condition), path("${condition}_symlink")
      tuple(val(condition), file(optim_var))
      tuple(val(condition), file(optim_q_qtl_bin))
    output:
      path("${outpath}/trans-by-cis_bonf_fdr.tsv", emit: trans_res, optional: true)
    script:

      // Use dosage?
      if (params.TensorQTL.use_gt_dosage) {
        dosage = "--dosage"
      }else{
        dosage = ""
      }

      // Define alpha value
      alpha = "0.05"

      // Generate the sumstats path
      sumstats_path = "${params.outdir}/TensorQTL_eQTLS/${condition}/"
        if (params.TensorQTL.interaction_file?.trim()) {
          inter_name = file(params.TensorQTL.interaction_file).baseName
          outpath_end = "interaction_output__${inter_name}"
      } else {
          inter_name = "NA"
          outpath_end = "base_output__base"
          }
      outpath = "./OPTIM_pcs/${outpath_end}"
    
      // Get the top result
      def value = ${optim_var}.text.trim()
      optim_value = value + "pcs"

      // Define the qval results file
      // qval_file="${sumstats_path}${nr_phenotype_pcs}/base_output/base/Cis_eqtls_qval.tsv"

      // Execute script, defininig whether or not it runs within the shell commands
      """
      echo ${optim_value} > temp.txt
      echo ${nr_phenotype_pcs} > temp2.txt
      if [ "${optim_value}" = "${nr_phenotype_pcs}" ]; then
        echo "yes" >> trans_run.txt
      else
        echo "not run, not optimum" >> trans_run.txt
      fi
      """

      // Execute if top result nPC and use value are the same    
      if (nr_phenotype_pcs == optim_value) {
        // Specify whether trans is ran
        """
        echo ${nr_phenotype_pcs} >> trans_run.txt
        echo "yes" >> trans_run.txt
        """
        // Run trans test
        """
        tensor_analyse_trans_by_cis.py \
          --covariates_file ${covariates_tsv} \
          --expression_bed Expression_Data.bed.gz \
          --plink_prefix_path ${plink_files_prefix}/plink_genotypes \
          --outdir ${outpath} \
          --dosage ${dosage} \
          --maf "0.05" \
          --cis_qval_results ${optim_q_qtl_bin} \
          --alpha ${alpha} \
          --window ${params.windowSize}
        
        cp -r ${condition}_symlink/"\$var"pcs__${outpath_end}/* ${outpath}
        """
      } else {
        """
        echo "not run, not optimum" >> trans_run.txt
        """
      }
}
*/

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
          if(params.TensorQTL.trans_by_cis){
            log.info 'Running trans-by-cis analysis on optimum nPCs'
            TRANS_BY_CIS(
              condition_bed,
              plink_genotype,
              OPTIMISE_PCS.out.optim_q_qtl_bin,
              OPTIMISE_PCS.out.optim_covariates,
              OPTIMISE_PCS.out.out_path
            )
          } 
  }
}
