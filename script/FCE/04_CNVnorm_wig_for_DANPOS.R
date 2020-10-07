args<-commandArgs(TRUE)

library(rtracklayer)

working.path<-args[1]
setwd(working.path)

out.path<-"analysis/03_CNV/wig_CNVnorm/"
dir.create(out.path)

#freec data
freec<-read.delim(grep("CNV", list.files("/analysis/03_CNV/"), value=T),header=F)
colnames(freec)<-c("chr", "start", "end", "CNV", "type")
freec[,1]<-paste("chr", freec[,1], sep="")

normWig<-function(gr.data, freec.fkt=freec){
        for(i in 1:nrow(freec.fkt)){
                chr<-freec.fkt$chr[i]
                start.iv<-(freec.fkt$start[i]+10)/10
                end.iv<-(freec.fkt$end[i]+1)/10
                len<-length(gr.data[[chr]])
                if(end.iv>len){end.iv=len}
                adj.score<-(gr.data[[chr]][start.iv:end.iv]$score)*(3/freec.fkt$CNV[i])
                gr.data[[chr]][start.iv:end.iv]$score<-adj.score
                print( paste(i,"out of",nrow(freec.fkt)) )
        }
        return(gr.data)
}


#load nuc occ data
file.path<-"analysis/02_preprocessing/bedpe/nucs/DANPOS/"

#get wig file
name.wig<-grep(".wig",list.files(paste0(file.path,"high_1n_r1/pooled/")), value=T)
wig.high.1r<-import.wig(paste0(file.path,"high_1n_r1/pooled/",name.wig),genome="hg19")

gr.list.high.1r<-split(wig.high.1r, seqnames(wig.high.1r))

gr.list.high.1r.norm<-normWig(gr.list.high.1r)
unlist.high.1r.norm<-unlist(gr.list.high.1r.norm)
export.wig(unlist.high.1r.norm, con=paste0(out.path,"high.1n.1r.norm.wig"))


#### low-MNase sample
name.wig<-grep(".wig",list.files(paste0(file.path,"low_1n_r1/pooled/")), value=T)
wig.low.1r<-import.wig(paste0(file.path,"low_1n_r1/pooled/",name.wig),genome="hg19")

gr.list.low.1r<-split(wig.high.1r, seqnames(wig.low.1r))
gr.list.low.1r.norm<-normWig(gr.list.low.1r)
unlist.low.1r.norm<-unlist(gr.list.low.1r.norm)
export.wig(unlist.low.1r.norm, con=paste0(out.path,"low.1n.1r.norm.wig"))




