args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-args[2]

## working path
setwd(working.path)
load("domain.df.rda")

dir.create("norm_GC/loess_norm/" ,recursive = T)


#libraries
library(ggplot2)
library(rtracklayer)
library("GenomicFeatures")
library(RColorBrewer)
library(colorRamps)
library(LSD)

#get colors
spec<-brewer.pal(10,"Spectral")

#filter mappability and extrem GC content
mapp.idx<-(domain.df$mappability> 0.9)
gc.idx<-domain.df$GC_content > 0.3 & domain.df$GC_content< 0.7

idx<-mapp.idx & gc.idx

### GC content
gc.domains<-domain.df$GC_content[idx]

##calculate ratio
high.r1<-domain.df$high_r1[idx]
low.r1<-domain.df$low_r1[idx]

fc.r1.high_low<-high.r1/low.r1


## filter
low.zero<-which(low.r1==0)
high.zero<-which(high.r1==0)
zeros<-unique(c(low.zero, high.zero))

fc.r1.high_low.cor<-high.r1[-zeros]/low.r1[-zeros]
gc.domains.cor<-gc.domains[-zeros]

ord<-order(gc.domains.cor)


####### loess normalization
high.loess <- loess(y ~ x, span=0.025, data.frame(x=gc.domains.cor, y=high.cor),
                       control=loess.control(surface = "interpolate", statistics = "none",
                                             cell=0.1))
high.predict <- predict(high.loess, data.frame(x=gc.domains.cor))


png("norm_GC/loess_norm/raw_withLoess_high.png")
    heatscatter(x=gc.domains.cor, y=high.cor,colpal=rev(spec), cor=T,
                cexplot=1, ylab="log2 high", xlab="GC content")
    lines((gc.domains.cor)[ord],high.predict[ord], col="black", lwd=2)
dev.off()

png("norm_GC/loess_norm/rawLOG_withLoess_high.png")
    heatscatter(x=gc.domains.cor, y=log2(high.cor),colpal=rev(spec), cor=T,
                cexplot=1, ylab="log2 high", xlab="GC content")
    lines((gc.domains.cor)[ord],log2(high.predict[ord]), col="black", lwd=2)
dev.off()


png("norm_GC/loess_norm/LoessNormed_high.png")
    heatscatter(x=gc.domains.cor, y=(high.cor)/high.predict,colpal=rev(spec), cor=T,
                cexplot=1, ylab="loessNorm log2 high", xlab="GC content")
    abline(h=1)
dev.off()

save(high.loess, file="norm_GC/loess_norm/high.loess.rda")




## low

low.loess <- loess(y ~ x, span=0.025, data.frame(x=gc.domains.cor, y=low.cor),
                    control=loess.control(surface = "interpolate", statistics = "none",
                                          cell=0.1))
low.predict <- predict(low.loess, data.frame(x=gc.domains.cor))


png("norm_GC/loess_norm/raw_withLoess_low.png")
    heatscatter(x=gc.domains.cor, y=low.cor,colpal=rev(spec), cor=T,
                cexplot=1, ylab="log2 low", xlab="GC content")
    lines((gc.domains.cor)[ord],low.predict[ord], col="black", lwd=2)
dev.off()

png("norm_GC/loess_norm/rawLOG_withLoess_low.png")
    heatscatter(x=gc.domains.cor, y=log2(low.cor),colpal=rev(spec), cor=T,
                cexplot=1, ylab="log2 low", xlab="GC content")
    lines((gc.domains.cor)[ord],log2(low.predict[ord]), col="black", lwd=2)
dev.off()


png("norm_GC/loess_norm/LoessNormed_low.png")
    heatscatter(x=gc.domains.cor, y=(low.cor)/low.predict,colpal=rev(spec), cor=T,
                cexplot=1, ylab="loessNorm log2 low", xlab="GC content")
    abline(h=1)
dev.off()

save(low.loess, file="norm_GC/loess_norm/low.loess.rda")


png("norm_GC/loess_norm/LoessNormed_logFC.png")
    heatscatter(y=log2((high.cor/high.predict)/(low.cor/low.predict)),
                x=gc.domains.cor,colpal=rev(spec), cor=T,
                cexplot=1, ylab="loessNorm log2 high", xlab="GC content")
    abline(h=0)
dev.off()



########## get Normalization factors from LOESS ###########


predicted.high<-predict(high.loess,data.frame(x=domain.df$GC_content))
predicted.low<-predict(low.loess,data.frame(x=domain.df$GC_content))

#filter
predicted.high[is.na(predicted.high)]<-0
predicted.low[is.na(predicted.low)]<-0

predicted.high[idx]<-0
predicted.low[idx]<-0

domain.df.high<-data.frame(domain.df[,1:3], normFactor=predicted.high)
domain.df.low<-data.frame(domain.df[,1:3], normFactor=predicted.low)



write.table(domain.df.low, file="norm_GC/loess_norm/low_norm_table.txt",
            sep="\t", quote=F, row.names = F)
write.table(domain.df.high, file="norm_GC/loess_norm/high_norm_table.txt",
            sep="\t", quote=F,row.names = F)