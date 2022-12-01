#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hailey.napier@duke.edu

#Hailey Napier
#November 30, 2022

cd /work/hcn4/nlgn2_mintchip_wd

#Load Java
module load Java/11.0.8 

#Execute BinarizeBam
java -mx4000M -jar ChromHMM/ChromHMM.jar BinarizeBam ChromHMM/CHROMSIZES/mm10.txt bams binarization_files/concatenated/nlgn2_concatenated_bam_binarization.tsv chromhmm_output_20221130/binarized

#Model Learning
java -mx4000M -jar ChromHMM/ChromHMM.jar LearnModel chromhmm_output_20221130/binarized chromhmm_output_20221130/learned_model 15 mm10


