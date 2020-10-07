AnalysisDir=$1
AnnoDir=$2
DataDir=$3

OutDir=$AnalysisDir/analysis/03_CNV

mkdir -p $OutDir


## map sonicated chromatin samples
cd $DataDir

bowtie2 -p 12 --very-sensitive-local \
-x $AnnoDir"/Bowtie2Index/genome" \
-U "input/"*.fastq.gz \
| samtools view -Sb -h -@ 12 -q 20 - > "input/input.bam"

samtools sort -@ 12  -o "input/input_sorted.bam" "input/input.bam"

rm "input/input.bam"

### run FREEC
~/tools/freec -conf $AnalysisDir"/script/additionalFiles/FREEC_config"
