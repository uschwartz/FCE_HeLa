AnalysisDir=$1
AnnoDir=$2

OutDir=$AnalysisDir/analysis/02_preprocessing
DataDir=$AnalysisDir/analysis/01_alignment

mkdir -p $OutDir"/bedpe/nucs/DANPOS"

### sort and size selection of fragments
cd $DataDir

for name in *_*; do

  echo 'sort and size selection of sample' $name

  #sort by read name
  samtools sort  -@ 12 -n -o $OutDir"/bedpe/"$name"_sorted_n.bam" $name".bam"

  cd $OutDir"/bedpe"

  #convert to bed in paired end format
  bamToBed -bedpe -i $name"_sorted_n.bam" >$name"_bedpe.bed"

  #remove sorted bam file
  rm $name"_sorted_n.bam"

  #select fragment sizes 140-200 bp
  awk  '$6-$2 > 139 && $6-$2 < 201' $name"_bedpe.bed" > "nucs/"$name"_nucs.bed"

  cd nucs

  #convert to bam and back to normal bed format
  bedtools bedpetobam -i  $name"_nucs.bed" -g $AnnoDir/ChromInfo.txt \
   >$name"_conversion.bam"
  bedtools bamtobed -i $name"_conversion.bam" >"DANPOS/"$name"_DANPOS.bed"

  rm $name"_conversion.bam"

  # run DANPOS2 to get nucleosome profiles
  cd DANPOS

  mkdir -p $name

  python /Users/admin/danpos-2.2.2/danpos.py dpos $name"_DANPOS.bed" -m 1 \
   --extend 70 -u 1e-2 -o $name >$name"/running_info.txt"


  cd $DataDir

done
