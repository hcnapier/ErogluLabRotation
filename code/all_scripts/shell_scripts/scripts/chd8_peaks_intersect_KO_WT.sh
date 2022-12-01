#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hailey.napier@duke.edu

# Hailey Napier
# December 1, 2022

#Assign directory variables
chd8_in_dir="/work/hcn4/nlgn2_mintchip_wd/chd8_peaks"
ko_wt_dir="/work/hcn4/nlgn2_mintchip_wd/chromhmm_output_20221130/learned_model"
out_dir="/work/hcn4/nlgn2_mintchip_wd/chd8_overlap/output_bed_files"
cd /work/hcn4/nlgn2_mintchip_wd/chd8_overlap

#Set up software
module load Bedtools/2.30.0  
homer="/work/hcn4/nlgn2_mintchip_wd/homer/bin"

#Convert .txt file of chd8 peaks to .bed file
perl $homer/pos2bed.pl $chd8_in_dir/merged_chd8_1_peaks_igg_input.txt_chd8_2_peaks_igg_input.txt > $chd8_in_dir/merged_chd8_peaks.bed

#Find intersection of chd8 peaks and KO segments
#KO input is segments .bed file because it has the start and end points as well as the state number, which is all we need to find overlap
#-a: file A
#-b: file B
#-u: Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B
bedtools intersect -a $ko_wt_dir/KO_15_segments.bed -b $chd8_in_dir/merged_chd8_peaks.bed -u > $out_dir/KO_overlap_chd8.bed

#Find intersection of chd8 peaks and WT segments
bedtools intersect -a $ko_wt_dir/WT_15_segments.bed -b $chd8_in_dir/merged_chd8_peaks.bed -u > $out_dir/WT_overlap_chd8.bed
