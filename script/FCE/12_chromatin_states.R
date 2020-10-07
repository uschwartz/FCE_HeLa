args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-args[2]
anno.path<-args[3]

## working path
setwd(working.path)

roadmap.out<-paste0(project.path,"/analysis/05_chromatin_features/chromStates")
dir.create(roadmap.out,recursive=T)

library(GenomicRanges)

#chrom length
chromInfo<-read.delim(file=paste0(anno.path,"/ChromInfo.txt"), header=F)
colnames(chromInfo)<-c("chr", "length")

##mappability
bigwig.rle<-import.bw(con=paste0(anno.path,
                "/mappability_hg19/wgEncodeCrgMapabilityAlign100mer.bigWig"),
                      as="RleList")

#chromatin states
segs<-read.delim(paste0(anno.path,"/E117_15_coreMarks_mnemonics.bed"), 
                 header=F)
colnames(segs)<-c("chr", "start", "end", "feature")

# subset
genome.seg<-segs[which(segs$chr %in% chromInfo$chr),]

seg.list<-split(genome.seg,f = genome.seg$feature)

names(seg.list)<-c("active_TSS","BivalentPoised_TSS","Flanking_Bivalent_EnhTSS","Bivalent_Enhancer",
                   "Repressed_PolyComb","Weak_Repressed_PolyComb", "Quiescent", "Flanking_Active_TSS",
                   "Transcribed_state_at_5_and_3","Strong_Transcription", "Weak_Transcription",
                   "Genic_Enhancer", "Enhancers", "ZNF_genes_and_repeats", "Heterochromatin")


save(seg.list, file=paste0(roadmap.out,"/seg.list.rda"))

############## get coverage at segments ######################

high.rle<-import.bw(con="high.1n.1r_GCnorm.bw",as="RleList")
low.rle<-import.bw(con="low.1n.1r_GCnorm.bw",as="RleList")



segments.anno.list<-list()
for(k in names(seg.list)){
    print(k)
    
    feature<-k  
    seg.feat<-seg.list[[feature]]
    chr<-as.character(unique(seg.feat$chr))
    
    
    feat.list<-list()
    for(i in chr){
        high.chr<-high.rle[[i]]
        low.chr<-low.rle[[i]]
        mapp.chr<-mapp.rle[[i]]
        
        chr.seg<-subset(seg.feat,(seg.feat$chr==i))
        agg.high<-aggregate(high.chr, start=chr.seg$start+1, end=chr.seg$end, FUN=mean)
        agg.low<-aggregate(low.chr, start=chr.seg$start+1, end=chr.seg$end, FUN=mean)
        
        agg.mapp<-aggregate(mapp.chr, start=chr.seg$start+1, end=chr.seg$end, FUN=mean)
        
        ##gc content
        print("read fasta")
        fasta<-read.fasta(file=paste0(anno.path,"/fasta/", chr.names, ".fa", sep=""), 
                          seqtype="DNA")
        sequence<-fasta[[1]] 
        print("fasta loaded")
        
        gc_count<-c()
        for(j in 1:nrow(chr.seg)){
            gc<-mean(sequence[chr.seg$start[j]:chr.seg$end[j]] %in% c("g", "c"))
            gc_count<-c(gc_count, gc)
        }
        
        feat.list[[i]]<-data.frame(chr.seg, high=agg.high,
                                   low=agg.low, mappability=agg.mapp,
                                   GC_content=gc_count)
        print(i)
    }
    
    feat.df<-do.call(rbind.data.frame, feat.list)

    segments.anno.list[[k]]<-feat.df
}

save(segments.anno.list, file=paste0(roadmap.out,"/segments.anno.list.rda"))


##################### plot data
features<-names(segments.anno.list)[c(10:15,1:9)]

col.seg<-c("firebrick","indianred", "salmon", "darkseagreen", "grey40","grey","white",
           "firebrick1", "springgreen", "springgreen4", "darkgreen", "darkorange","gold",
           "turquoise","slateblue")



feat.list<-list()
feat.list.gc<-list()
for( k in features){
    feat<-segments.anno.list[[k]]
    #filter mappability < 0.9
    feat.filt<-feat[feat$mappability>0.9,]
    #filter zeros
    

    zero.gc<-which(feat.filt$low==0 | feat.filt$high==0)
    
    
    if(length(zero) !=0){
        feat.list[[k]]<-data.frame(cond=rep(k, (nrow(feat.filt)-length(zero))),
                                   logFC=log2(feat.filt$high[-zero]/feat.filt$low[-zero]),
                                   GC=feat.filt$GC_content[-zero],
                                   width=feat.filt$end[-zero]-feat.filt$start[-zero])
    } else {
        feat.list[[k]]<-data.frame(cond=rep(k, (nrow(feat.filt))),
                                   logFC=log2(feat.filt$high/feat.filt$low),
                                   GC=feat.filt$GC_content,
                                   width=feat.filt$end-feat.filt$start)
    }
    
    feat.list.gc[[k]]<-data.frame(cond=rep(k, (nrow(feat.filt)-length(zero.gc))),
                                  logFC=log2(feat.filt$high[-zero.gc]/feat.filt$low[-zero.gc]))
    print(k)
}


feat.gc.df<-do.call(rbind.data.frame, feat.list.gc)
feat.df<-do.call(rbind.data.frame, feat.list)


library(ggplot2)

p<-ggplot(feat.gc.df, aes(cond,logFC))
p<-p+geom_boxplot(aes(fill=cond), outlier.shape = NA)+scale_y_continuous(limits=c(-2,2))+
    scale_fill_manual(values=col.seg)+geom_hline(yintercept = 0, linetype="dashed")

pdf(paste0(roadmap.out,"/boxplot_chromatin_states_logfc.pdf"), width=20)
    print(p)
dev.off()

## Width

p<-ggplot(feat.df, aes(cond,width))
p<-p+geom_boxplot(aes(fill=cond), outlier.shape = NA)+scale_fill_manual(values=col.seg)+
    scale_y_continuous(limits=c(0,10000))

pdf(paste0(roadmap.out,"/boxplot_chromatin_states_width.pdf"), width=20)
    print(p)
dev.off()
