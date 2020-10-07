args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-args[2]
anno.path<-args[3]
setwd(working.path)

#Libaries
library(rtracklayer)
library("GenomicFeatures")
library(seqinr)

#load data
chrom.sizes<-read.delim(file=paste0(anno.path,"/ChromInfo.txt"), header=F)

load("high.1r.rda")
load("low_1r.rda")

domain.list<-list()

for(i in chrom.sizes[,1]){
        domain.list[[i]]<-data.frame(chrom=high.1r[[i]]$chrom,
                                     start=high.1r[[i]]$start,
                                     end=high.1r[[i]]$end,
                                     high.1r=high.1r[[i]]$mean,
                                     low.1r=low.1r[[i]]$mean)
}

bigwig.rle<-import.bw(con=paste0(anno.path,
                     "/mappability_hg19/wgEncodeCrgMapabilityAlign100mer.bigWig"),
                      as="RleList")
GCcontent.list<-list()
mappab.list<-list()

for ( i in 1:nrow(chrom.sizes)){

    chr.names=as.character(chrom.sizes[i,1]) 
    chr.length<-chrom.sizes[i,2]
    domains.chr<-domain.list[[chr.names]]
    
    ##gc content
    print("read fasta")
    fasta<-read.fasta(file=paste0(anno.path,"/fasta/", chr.names, ".fa", sep=""), seqtype="DNA")
    sequence<-fasta[[1]] 
    print("fasta loaded")
    
    gc_count<-c()
    for(j in 1:nrow(domains.chr)){
        gc<-mean(sequence[domains.chr$start[j]:domains.chr$end[j]] %in% c("g", "c"))
        gc_count<-c(gc_count, gc)
   }
    
    mappability<-c()
    bigwig.chr<-bigwig.rle[[chr.names]]
    
    slid.win.mean<-aggregate(bigwig.chr, start=domains.chr$start,
                             end=domains.chr$end, FUN=mean)
    valid=slid.win.mean>0.9
    
    domain.list[[chr.names]]<-data.frame(domain.list[[chr.names]], GC_content=gc_count,
                                         mappability=slid.win.mean, valid)
    GCcontent.list[[chr.names]]<-data.frame(domains.chr[,1:3],GC_content=gc_count)
    mappab.list[[chr.names]]<-data.frame(domains.chr[,1:3],mappability=slid.win.mean)
    
    print(chr.names)
}



domain.df<-do.call(rbind.data.frame, domain.list)

median.sig<-apply(domain.df[domain.df$valid,c(4:7)],1, median)
sd<-sd(median.sig)

pdf("filtering_3sd.pdf")
    plot(density(median.sig), xlab="median signal", main="filter deletions and dulications")
    abline(v=mean(median.sig)-3*sd)
    abline(v=mean(median.sig)+3*sd)
dev.off()


invalid<-(median.sig< mean(median.sig)-3*sd)  | (median.sig > mean(median.sig)+3*sd) 

val.vec<-domain.df$valid
val.vec[val.vec][invalid]<-rep(FALSE,sum(invalid))
domain.df$valid<-val.vec

save(domain.df, file="domain.df.rda")

domain.list<-split(domain.df,domain.df$chrom)

save(GCcontent.list, file="GCcontent.list.rda")
save(mappab.list, file="mappab.list.rda")
save(domain.list, file="domain.list.rda")




