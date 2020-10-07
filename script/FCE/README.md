## Content

+ 00_execute_scripts.sh : Masterscript to configure project and annotation path and executes the following scripts
+ 01_align2genome.sh : Bash-script to map raw MNase-seq data to the human genome
+ 02_bam2nucleosomes.sh : Bash-script to process properly aligned read pairs and to generate nucleosome frequency maps
+ 03_CNV_correction.sh : Bash-script to call CNV in HeLa cells; used in combination with "../additionalFiles/FREEC_config"
