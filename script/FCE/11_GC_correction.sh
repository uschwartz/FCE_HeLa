AnalysisDir=$1
AnnoDir=$2
DataDirGC=$3

DataDirWig=$AnalysisDir"/analysis/03_CNV/wig_CNVnorm/"


cd $DataDirGC

$AnalysisDir"/script/additionalFiles//GC-Parser/Parser" $DataDirWig"/high.1n.1r.norm.wig" high_norm_table.txt high.1n.1r_GCnorm.wig

$AnalysisDir"/script/additionalFiles//GC-Parser/Parser" $DataDirWig"/low.1n.1r.norm.wig" low_norm_table.txt low.1n.1r_GCnorm.wig

# wig2bigwig

wigToBigWig high.1n.1r_GCnorm.wig \
$AnnoDir"/ChromInfo.txt" high.1n.1r_GCnorm.bw

wigToBigWig low.1n.1r_GCnorm.wig \
 $AnnoDir"/ChromInfo.txt" low.1n.1r_GCnorm.bw
