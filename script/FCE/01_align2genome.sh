AnalysisDir=$1
AnnoDir=$2
DataDir=$3

OutDir=$AnalysisDir/analysis/01_alignment
OutDirQC=$AnalysisDir/QC/01_alignment


### alignment
mkdir -p $OutDir
mkdir -p $OutDirQC


cd $DataDir

for entry in *_*; do

  echo 'align to genome' $entry

  mkdir $OutDir"/"$entry

	bowtie2 -p 12 --very-sensitive-local --no-discordant \
  -x $AnnoDir"/Bowtie2Index/genome" \
  -1 $entry"/"*pair1.fastq.gz -2 $entry"/"*pair2.fastq.gz \
  2> $OutDirQC'/'$entry"-align-log.txt" \
  | samtools view -Sb -h -@ 12 -f 2 -q 20 - > $OutDir"/"$entry"/"$entry".bam"

done
