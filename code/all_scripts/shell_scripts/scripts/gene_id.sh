#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hailey.napier@duke.edu

# Hailey Napier
# December 13, 2022

#Assign directory variables
state_overlap_in="/work/hcn4/nlgn2_mintchip_wd/gene_id/state_overlap_files"
mm10_known_genes="/work/hcn4/nlgn2_mintchip_wd/gene_id/mm10.knownGene.gtf"
mm10_refseq="/work/hcn4/nlgn2_mintchip_wd/gene_id/mm10.ncbiRefSeq.gtf" 
out_dir="/work/hcn4/nlgn2_mintchip_wd/gene_id/overlap_out/refseq"
cd /work/hcn4/nlgn2_mintchip_wd/gene_id/state_overlap_files

#Set up software
module load Bedtools/2.30.0  

#Find intersection of chromhmm chd8 binding sites and gene id info
#KO E4
#-a: file A
#-b: file B
#-u: Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B
#bedtools intersect -b $state_overlap_in/overlap_KO_E4.bed -a $mm10_known_genes -u > $out_dir/KO_E4_known_genes_chromhmm_chd8_bind.bed
#KO E5
#bedtools intersect -b $state_overlap_in/overlap_KO_E5.bed -a $mm10_known_genes -u > $out_dir/KO_E5_known_genes_chromhmm_chd8_bind.bed
#WT E4
#bedtools intersect -b $state_overlap_in/overlap_WT_E4.bed -a $mm10_known_genes -u > $out_dir/WT_E4_known_genes_chromhmm_chd8_bind.bed
#WT E5
#bedtools intersect -b $state_overlap_in/overlap_WT_E5.bed -a $mm10_known_genes -u > $out_dir/WT_E5_known_genes_chromhmm_chd8_bind.bed

# Try again with refseq
#KO E4
bedtools intersect -b $state_overlap_in/overlap_KO_E4.bed -a $mm10_refseq -u > $out_dir/KO_E4_refseq_chromhmm_chd8_bind.bed
#KO E5
bedtools intersect -b $state_overlap_in/overlap_KO_E5.bed -a $mm10_refseq -u > $out_dir/KO_E5_refseq_chromhmm_chd8_bind.bed
#WT E4
bedtools intersect -b $state_overlap_in/overlap_WT_E4.bed -a $mm10_refseq -u > $out_dir/WT_E4_refseq_chromhmm_chd8_bind.bed
#WT E5
bedtools intersect -b $state_overlap_in/overlap_WT_E5.bed -a $mm10_refseq -u > $out_dir/WT_E5_refseq_chromhmm_chd8_bind.bed
