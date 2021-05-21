#!/bin/bash
#SBATCH --job-name=motevo
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --output=./log_mot/%A_%a.out
#SBATCH --error=./log_mot/%A_%a.err
#SBATCH --mem=4G
#SBATCH --array=1-110
      
mkdir -p log
mkdir -p log_mot

perl /scicore/home/nimwegen/GROUP/AlignmentPipeline/scripts/motevoc_wrapper.pl --alignments ./list_of_prom.fasta --wm_dir /scicore/home/nimwegen/urchuegu/projects/e_coli_mara/annotation/weight_matrices/rdb_updated/individual_wms --working_dir . --motevo_path /scicore/home/nimwegen/GROUP/software/motevo_ver1.11/bin/motevo  &> log/motevo.$SLURM_ARRAY_TASK_ID

 
