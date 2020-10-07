AnalysisDir=$1


DataDir<-$AnalysisDir"/analysis/03_CNV/wig_CNVnorm/"
OutDir<-$AnalysisDir"/analysis/03_CNV/wig_CNVnorm/meta_profiles"
mkdir -p $OutDir

## DNase
mkdir $OutDir"/DNase"
cd $OutDir"/DNase"

python /Users/admin/danpos-2.2.2/danpos.py profile \
 $DataDir"/high_1n_r1/pooledadjust_Karyotype_high.1n.1r.norm.smooth.wig",\
$DataDir"/low_1n_r1/pooledadjust_Karyotype_low.1n.1r.norm.smooth.wig" \
 --bed3file_paths $AnalysisDir"/script/additionalFiles/meta-plots/ENCFF692NCU_DNase.bed" \
  --flank_up 2000  --flank_dn 2000 --genomic_sites center \
   --plot_column 1 --plot_row 1 --wigfile_aliases high,low

## enhancer
mkdir $OutDir"/enhancer"
cd $OutDir"/enhancer"

python /Users/admin/danpos-2.2.2/danpos.py profile \
  $DataDir"/high_1n_r1/pooledadjust_Karyotype_high.1n.1r.norm.smooth.wig",\
  $DataDir"/low_1n_r1/pooledadjust_Karyotype_low.1n.1r.norm.smooth.wig" \
  --bed3file_paths $AnalysisDir"/script/additionalFiles/meta-plots/enhancer_anderson.bed" \
  --flank_up 2000  --flank_dn 2000 --genomic_sites center \
  --plot_column 1 --plot_row 1 --wigfile_aliases high,low



## genes
mkdir $OutDir"/genes"
cd $OutDir"/genes"

python /Users/admin/danpos-2.2.2/danpos.py profile \
  $DataDir"/high_1n_r1/pooledadjust_Karyotype_high.1n.1r.norm.smooth.wig",\
  $DataDir"/low_1n_r1/pooledadjust_Karyotype_low.1n.1r.norm.smooth.wig" \
  --genefile_paths $AnalysisDir"/script/additionalFiles/meta-plots/ENCFF000DNW_highly_expressed.txt" \
  --flank_up 2000  --flank_dn 2000 --genomic_sites TSS,gene \
  --plot_column 1 --plot_row 1 --wigfile_aliases high,low
