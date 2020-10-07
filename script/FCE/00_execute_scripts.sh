######## specify path to project folder an annotation infromation ##############

#project folder
AnalysisDir=~/Analysis/FCE_HeLa/

#folder with bowtie2 index
AnnoDir=~/Analysis/FCE_HeLa/data/annotation

#folder with raw data
rawDir=~/Analysis/FCE_HeLa/data/sequencing

################################################################################
# Alignment
$AnalysisDir/script/FCE/01_align2genome.sh $AnalysisDir $AnnoDir $rawDir

# get Nucleosome profiles
$AnalysisDir/script/FCE/02_bam2nucleosomes.sh $AnalysisDir $AnnoDir

# call CNV
$AnalysisDir/script/FCE/03_CNV_correction.sh $AnalysisDir $AnnoDir $rawDir

# normalize
Rscript $AnalysisDir/script/FCE/04_CNVnorm_wig_for_DANPOS.R $AnalysisDir

# get coverage values  in bins
Rscript $AnalysisDir/script/FCE/05_window_analysis.R $AnalysisDir

# get GC and mappability in bins
Rscript $AnalysisDir/script/FCE/06_get_mappability_and_GC.R $AnalysisDir \
 $AnalysisDir"analysis/04_bins/win1e+06" $AnnoDir

Rscript $AnalysisDir/script/FCE/06_get_mappability_and_GC.R $AnalysisDir \
  $AnalysisDir"analysis/04_bins/win10000" $AnnoDir

Rscript $AnalysisDir/script/FCE/06_get_mappability_and_GC.R $AnalysisDir \
   $AnalysisDir"analysis/04_bins/win1000" $AnnoDir

Rscript $AnalysisDir/script/FCE/06_get_mappability_and_GC.R $AnalysisDir \
    $AnalysisDir"analysis/04_bins/win500" $AnnoDir

Rscript $AnalysisDir/script/FCE/06_get_mappability_and_GC.R $AnalysisDir \
        $AnalysisDir"analysis/04_bins/win250" $AnnoDir

### correlate GC content and nucleosome frequency
Rscript $AnalysisDir/script/FCE/07_plot_GC_correlation.R $AnalysisDir \
        $AnalysisDir"analysis/04_bins/win1e+06"

Rscript $AnalysisDir/script/FCE/07_plot_GC_correlation.R $AnalysisDir \
            $AnalysisDir"analysis/04_bins/win250"


### call nucleosome positions
$AnalysisDir/script/FCE/08_nucleosome_calling.sh $AnalysisDir $AnnoDir

#nucleosome meta-profiles at DNase HS; enhancer; fit_regression
$AnalysisDir/script/FCE/09_meta_profiles.sh $AnalysisDir

### get loess function for GC normlization
Rscript $AnalysisDir/script/FCE/10_get_loess_GC_fkt.R $AnalysisDir \
            $AnalysisDir"analysis/04_bins/win250"

# run GC normalization
$AnalysisDir/script/FCE/11_GC_correction.sh $AnalysisDir $AnnoDir \
  $AnalysisDir"analysis/04_bins/win250/norm_GC/loess_norm/"

## chromatin states
Rscript $AnalysisDir/script/FCE/12_chromatin_states.R $AnalysisDir \
  $AnalysisDir"analysis/04_bins/win250/norm_GC/loess_norm/" $AnnoDir


## histone modifications
Rscript $AnalysisDir/script/FCE/13_histone_marks.R $AnalysisDir \
    $AnalysisDir"analysis/04_bins/win250/norm_GC/loess_norm/" $AnnoDir
