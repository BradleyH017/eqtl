#! bin/bash


# Remove old logs but not the most previous run
# rm -r *html.*;
# rm .nextflow.log.*;
# rm flowchart.png.*;
# rm trace.txt.*;

# TODO: User edit this section to fit workstation
# export NF_PATH="/software/hgi/installs/anaconda3/envs/nextflow/bin/"
export REPOS="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/"
export EQTL_REPO="${REPOS}/eqtl"
export INPUT="${EQTL_REPO}/sample_input/TI_fr003_v004-input/input_minimal.nf"

# Add Singularity to path
# PATH=$PATH:/software/singularity/v3.10.0/bin

module load ISG/singularity/3.6.4
module load HGI/common/nextflow/22.04.4
singularity --version

# Nextflow settings
# export NXF_OPTS="-Xms25G -Xmx25G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
# export NXF_HOME=$(pwd)
# export NXF_WORK="${NXF_HOME}/work"
# export NXF_TEMP="${NXF_HOME}/nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/nextflow_conda"
export SINGULARITY_CACHEDIR="${NXF_HOME}/cache_singularity"
# export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/cache_singularity"
# Farm specific settings
# export QT_QPA_PLATFORM='offscreen'
# export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# export LSB_DEFAULTGROUP="sc-eqtl-ibd"
# export LSB_DEFAULTGROUP="cnv_15x"
# export LSB_DEFAULTGROUP="team152"
# export PYTHONHASHSEED=0

# # Uncomment this if get strange bus errors
# # export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality

# # mkdir $NXF_TEMP
# # mkdir results
# # cp ${INPUT} results/
# # cp ${EQTL_REPO}/conf/extra_confs/* results/

# # To stop TclError
# # export QT_QPA_PLATFORM='offscreen'

# Run nextflow
nextflow run ${EQTL_REPO} -profile sanger -c ${INPUT} -resume

# Execute:
# module load ISG/singularity/3.6.4
# module load HGI/common/nextflow/22.04.4
# sh run_nfeqtl.sh