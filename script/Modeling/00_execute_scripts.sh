######## specify path to project folder an annotation infromation ##############

#project folder
AnalysisDir=~/Analysis/FCE_HeLa/


################################################################################
# simulate MNase ladder
Rscript $AnalysisDir/script/Modeling/01_getLadder.R $AnalysisDir

#simulate output of MNase digestion
Rscript $AnalysisDir/script/Modeling/02_simulate_FCE.R $AnalysisDir
