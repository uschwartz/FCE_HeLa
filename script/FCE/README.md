## Content

+ 00_execute_scripts.sh : Masterscript to configure project and annotation path and executes the following scripts
+ 01_align2genome.sh : Bash-script to map raw MNase-seq data to the human genome
+ 02_bam2nucleosomes.sh : Bash-script to process properly aligned read pairs and to generate nucleosome frequency maps
+ 03_CNV_correction.sh : Bash-script to call CNV in HeLa cells; used in combination with config file "../additionalFiles/FREEC_config"
+ 04_CNVnorm_wig_for_DANPOS.R : R script to normalize nucleosome frequency maps by CNVs and generates a wig output file
+ 05_window_analysis.R : R script to segement the genome into bins and assess the normalized nucleosome abundence of each sample within each bin
+ 06_get_mappability_and_GC.R : R script to assess average GC content and mappability score within bins
+ 07_plot_GC_correlation.R : R script to assess the GC dependency of each MNase digestion condition 
+ 08_nucleosome_calling.sh : Bash-script to call nucleosome positions onto the CNV normalized nucleosome frequency maps
+ 09_meta_profiles.sh : Bash-script to generate meta profiles of nucleosome occupancy around: enhancer, DNase HS sites, TSS of transcribed genes
+ 10_get_loess_GC_fkt.R : R script to calculate the GC normalization factors using loess curve fitting
+ 11_GC_correction.sh : Bash-script to normalize the MNase GC dependency in the nucleosome frequency profile; used in combination with in house script "../additionalFiles/GC-parser/Parser"
+ 12_chromatin_states.R : R script to assess the GC normalized nucleosome abundance in annotated chromatin states
+ 13_histone_marks.R : R script to correlate histone ChIP-seq read density with GC normalized nucleosome abundance
