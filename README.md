# ErogluLabRotation
Hailey Napier
Duke University
10/22/22 - 12/16/22

This repository will serve as documentation for the computational aspects of my rotation project in the Eroglu Lab (10/22/22 - 12/16/22). 

See my labarchives folder ("Hailey N. Oct 2022" in the "Rotations Students" directory within Kristina Sakers's labarchives notebook) for full descriptions of experiments and analysis. 

### Directory and File Descriptions
+ **r_scripts** contains only r scripts
  + **chromhmm_scripts** contains r scripts used to analyze and create plots for chromhmm data and data related to chromatin states/chd8 binding 
    + *chd8_overlap_wt_ko.R* is a script to find overlap between chd8 binding sites and chromatin states
    + *get_chd8_binding_genes.R* is a script to find the genes chd8 binds and their associated GO terms
    + *make_bam_binarization_files.R* is a script to format the files for binarization for chromhmm
  + **image_analysis** contains r scripts used to quantify image analysis and make plots for image analysis data
    + *20221206_Cx43_image_analysis.R* is a script to analyze connexin43 puncta
+ **shell_scripts** contains shell scripts I ran on DCC
  + *chd8_peaks_intersect_KO_WT.sh* is a script to find the overlap of chd8 binding with chromatin states
  + *nlgn2_chromhmm.sh* is a script to find chromatin states using mintchip .bed files from nlgn2 KO astrocytes
+ *backup_scripts.sh* is a script to backup DCC shell files on my local computer (to be run on local terminal, not DCC)


Contact hailey.napier@duke.edu with questions.
