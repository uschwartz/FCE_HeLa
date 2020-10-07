args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-args[2]
anno.path<-args[3]

## working path
setwd(working.path)
#output to 
out.path<-"histoneMarks"
dir.create(out.path)

# libraries
library(rtracklayer)
library("GenomicFeatures")
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(LSD)

#chrom length
chrom.sizes<-read.delim(file=paste0(anno.path,"/ChromInfo.txt"), header=F)
#mappability
bigwig.rle<-import.bw(con=paste0(anno.path,
                    "/mappability_hg19/wgEncodeCrgMapabilityAlign100mer.bigWig"),
                      as="RleList")

getNucOccupancy<-function(data, name.exp ,win=1e6,chromosomes=chrom.sizes){
    
    domains.list<-vector("list",nrow(chromosomes))
    names(domains.list)<-as.character(chromosomes[,1])
    
    
    for(i in 1:nrow(chromosomes)){
        #parameters
        chrom.name<-as.character(chromosomes[i,1])
        endpoint<-chromosomes[i,2]
        
        #select chromosome
        chr.gr<-data[[chrom.name]]
        
        #convert to Rle object
        chr.rle<-Rle(chr.gr$score, width(chr.gr))
        
        #fill up with zeros
        app<-rep(0,endpoint-length(chr.rle))
        chr.rle.ext<-c(chr.rle,app)
        
        #calculate average, median, variation and sum
        slid.win.mean<-aggregate(chr.rle.ext, start=seq(1,endpoint-win,win),
                                 width=win, FUN=mean)
        slid.win.median<-aggregate(chr.rle.ext, start=seq(1,endpoint-win,win),
                                   width=win, FUN=median)
        slid.win.sum<-aggregate(chr.rle.ext, start=seq(1,endpoint-win,win),
                                width=win, FUN=sum)
        slid.win.var<-aggregate(chr.rle.ext, start=seq(1,endpoint-win,win),
                                width=win, FUN=var)
        
        
        numb.win<-length(slid.win.mean)
        
        
        ##data frame
        df<-data.frame(chrom=rep(chrom.name,numb.win),start=seq(1,(numb.win-1)*win+1,win),
                       end=seq(win,numb.win*win,win), mean=slid.win.mean,
                       median=slid.win.median, sum=slid.win.sum, var=slid.win.var)
        
        domains.list[[chrom.name]]<-df
        print(chrom.name)
    }
    
    return(domains.list)
}


############## get coverage traccks ######################

high_normGC<-import.bw(con="high.1n.1r_GCnorm.bw", genome="hg19")
low_normGC<-import.bw(con="low.1n.1r_GCnorm.bw", genome="hg19")

#split chromosomes
gr.list.high_normGC<-split(high_normGC, seqnames(high_normGC))
gr.list.low_normGC<-split(low_normGC, seqnames(low_normGC))

varWindow<-function(win){
    dir.create(paste("win",win,sep=""))
    high.r1<-getNucOccupancy(data=gr.list.high_normGC, name.exp="high_r1",
                             win=win,chromosomes=chrom.sizes)
    save(high.r1,file=paste(out.path,"/high.r1.rda",sep=""))
    
    low.r1<-getNucOccupancy(data=gr.list.low_normGC, name.exp="low_r1",
                            win=win,chromosomes=chrom.sizes)
    save(low.r1,file=paste(out.path,"/low.r1.rda",sep=""))
}


##run analysis
print("start  win 1e4")
varWindow(1e4)
print("end win 1e4")

############## get ChIP inetnsities ##############
domain.list<-list()

for(i in chrom.sizes[,1]){
    domain.list[[i]]<-data.frame(chrom=high.r1[[i]]$chrom,
                                 start=high.r1[[i]]$start,
                                 end=high.r1[[i]]$end,
                                 high_r1=high.r1[[i]]$mean,
                                 low_r1=low.r1[[i]]$mean)
}


#load encode data
H3K27ac.rle<-import.bw(con=paste0(anno.path,
                        "/ENCFF311EWS_H3K27ac_Broad_FC_over_Ctrl.bigwig"),
                       as="RleList")
H3K27me3.rle<-import.bw(con=paste0(anno.path,
                        "/ENCFF958BAN_H3K27me3_Broad_FC_over_Ctrl.bigwig"),
                        as="RleList")


GCcontent.list<-list()
mappab.list<-list()
H3K27ac.list<-list()
H3K27me3.list<-list()

for ( i in 1:nrow(chrom.sizes)){
    
    chr.names=as.character(chrom.sizes[i,1]) 
    chr.length<-chrom.sizes[i,2]
    domains.chr<-domain.list[[chr.names]]
    
    ##gc content
    print("read fasta")
    fasta<-read.fasta(file=paste0(anno.path,"/fasta/", chr.names, ".fa", sep=""), 
                      seqtype="DNA")
    
    sequence<-fasta[[1]] 
    print("fasta loaded")
    
    gc_count<-c()
    for(j in 1:nrow(domains.chr)){
        gc<-mean(sequence[domains.chr$start[j]:domains.chr$end[j]] %in% c("g", "c"))
        gc_count<-c(gc_count, gc)
    }
    
    bigwig.chr<-bigwig.rle[[chr.names]]
    slid.win.mean<-aggregate(bigwig.chr, start=domains.chr$start,
                             end=domains.chr$end, FUN=mean)
    
    H3K27ac.chr<-H3K27ac.rle[[chr.names]]
    H3K27ac.mean<-aggregate(H3K27ac.chr, start=domains.chr$start,
                            end=domains.chr$end, FUN=mean)
    
    H3K27me3.chr<-H3K27me3.rle[[chr.names]]
    H3K27me3.mean<-aggregate(H3K27me3.chr, start=domains.chr$start,
                             end=domains.chr$end, FUN=mean)
    
    valid=slid.win.mean>0.9
    
    domain.list[[chr.names]]<-data.frame(domain.list[[chr.names]], GC_content=gc_count,
                                         mappability=slid.win.mean, valid,
                                         H3K27ac=H3K27ac.mean,H3K27me3=H3K27me3.mean )
    
    print(chr.names)
}



domain.df<-do.call(rbind.data.frame, domain.list)
domain.df.noCORR<-domain.df
save(domain.df.noCORR, file=paste0(out.path,"/domain.df.noCORR.rda"))


## plot
domain.df<-domain.df.noCORR[domain.df.noCORR$valid,c("high_r1", "low_r1", "H3K27ac", "H3K27me3")]

df.H3K27<-data.frame(logFC=log2(domain.df$high_r1/domain.df$low_r1),
                     H3K27ac=log2(domain.df$H3K27ac),
                     H3K27me3=log2(domain.df$H3K27me3)
    )

# get quantiles/quartiles
##################################### H3K27ac

wx.H3K27ac<-which(df.H3K27$H3K27ac>0)
q.H3K27ac<-quantile(df.H3K27$H3K27ac[wx.H3K27ac])

H3K27ac_quartiles<-rep("no Enrichment",nrow(df.H3K27))
H3K27ac_quartiles[df.H3K27$H3K27ac>0]<-"quartile1"
H3K27ac_quartiles[df.H3K27$H3K27ac>q.H3K27ac["25%"]]<-"quartile2"
H3K27ac_quartiles[df.H3K27$H3K27ac>q.H3K27ac["50%"]]<-"quartile3"
H3K27ac_quartiles[df.H3K27$H3K27ac>q.H3K27ac["75%"]]<-"quartile4"

##################################### H3K27me3

wx.H3K27me3<-which(df.H3K27$H3K27me3>0)
q.H3K27me3<-quantile(df.H3K27$H3K27me3[wx.H3K27me3])


H3K27me3_quartiles<-rep("no Enrichment",nrow(df.H3K27))
H3K27me3_quartiles[df.H3K27$H3K27me3>0]<-"quartile1"
H3K27me3_quartiles[df.H3K27$H3K27me3>q.H3K27me3["25%"]]<-"quartile2"
H3K27me3_quartiles[df.H3K27$H3K27me3>q.H3K27me3["50%"]]<-"quartile3"
H3K27me3_quartiles[df.H3K27$H3K27me3>q.H3K27me3["75%"]]<-"quartile4"


#### boxplot

df<-cbind(df.H3K27, H3K27me3_quartiles, H3K27ac_quartiles)

H3K27ac.part<-cbind(df[,c("logFC", "H3K27ac", "H3K27ac_quartiles")],rep("H3K27ac", nrow(df)))
colnames(H3K27ac.part)<-c("logFC", "HistoneEnrichment", "Quartiles", "Histone_modification")      

H3K27me3.part<-cbind(df[,c("logFC", "H3K27me3", "H3K27me3_quartiles")],rep("H3K27me3", nrow(df)))
colnames(H3K27me3.part)<-c("logFC", "HistoneEnrichment", "Quartiles", "Histone_modification")      

df.plot<-rbind(H3K27ac.part,H3K27me3.part)
df.plot$Histone_modification<-relevel(df.plot$Histone_modification, "H3K27me3")

g<-ggplot(df.plot, aes(x=Quartiles, y=logFC))
g<-g+geom_boxplot(aes(fill=Histone_modification), outlier.shape = NA)+
    scale_y_continuous(limits=c(-0.5,0.5))+
    scale_fill_manual(values=c("dodgerblue4","darkorange"))+geom_hline(yintercept = 0, linetype="dashed")

pdf(paste0(out.path,"/HistoneModifications_Quartiles.pdf"), width=10, height = 5)
    print(g)
dev.off()