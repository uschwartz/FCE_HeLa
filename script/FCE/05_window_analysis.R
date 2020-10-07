args<-commandArgs(TRUE)

library(rtracklayer)

working.path<-args[1]
setwd(working.path)

import.path<-"analysis/03_CNV/wig_CNVnorm/"
out.path<-"analysis/04_bins/"
dir.create(out.path)
####

chrom.sizes<-read.delim(file="data/annotation/ChromInfo.txt", header=F)

######

#load coverage profiles
#high-MNase
wig.high.1r<-import.wig(paste0(import.path,"high.1n.1r.norm.wig"),genome="hg19")
gr.list.high.1r<-split(wig.high.1r, seqnames(wig.high.1r))

#low-MNase
wig.low.1r<-import.wig(paste0(import.path,"low.1n.1r.norm.wig"),genome="hg19")
gr.list.low.1r<-split(wig.low.1r, seqnames(wig.low.1r))

#####

### functions

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


varWindow<-function(win){
    setwd(out.path)
    dir.create(paste("win",win,sep=""))
    #high-MNase
    high.1r<-getNucOccupancy(data=gr.list.high.1r, name.exp="high_1r",
                            win=win,chromosomes=chrom.sizes)
    save(high.1r,file=paste("win",win,"/high.1r.rda",sep=""))
    #low-MNase
    low.1r<-getNucOccupancy(data=gr.list.low.1r, name.exp="low_1r",
                             win=win,chromosomes=chrom.sizes)
    save(low.1r,file=paste("win",win,"/low.1r.rda",sep=""))

}


##run
print("start  win 1e4")
varWindow(1e4)
print("end win 1e4")

print("start  win 1e3 ")
varWindow(1e3)
print("end win 1e3")

print("start  win 500 ")
varWindow(500)
print("end win 500")

print("start  win 250 ")
varWindow(250)
print("end win 250")

print("start  win 1e6")
varWindow(1e6)
print("end win 1e6")
