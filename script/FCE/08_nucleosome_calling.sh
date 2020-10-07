AnalysisDir=$1
AnnoDir=$2

DataDir<-$AnalysisDir"/analysis/03_CNV/wig_CNVnorm/"


cd $DataDir

python /Users/admin/danpos-2.2.2/danpos.py dpos high.1n.1r.norm.wig   -o high_1n_r1

python /Users/admin/danpos-2.2.2/danpos.py dpos low.1n.1r.norm.wig   -o low_1n_r1


# wig2bigwig

cd $DataDir"/high_1n_r1/pooled"
wigToBigWig adjust_Karyotype_high.1n.1r.norm.smooth.wig \
 $AnnoDir"/ChromInfo.txt" adjust_Karyotype_high.1n.1r.norm.smooth.bw


cd $DataDir"/low_1n_r1/pooled"
wigToBigWig adjust_Karyotype_low.1n.1r.norm.smooth.wig \
$AnnoDir"/ChromInfo.txt" adjust_Karyotype_low.1n.1r.norm.smooth.bw
