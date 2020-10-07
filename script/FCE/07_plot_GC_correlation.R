args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-args[2]

## working path
setwd(working.path)
load("domain.list.rda")

## libraries
library(ggplot2)
library(rtracklayer)
library("GenomicFeatures")
library(RColorBrewer)
library(colorRamps)
library(LSD)


# colors
greens<-brewer.pal(9, "Greens")
reds<-brewer.pal(9, "Reds")
greys<-brewer.pal(9, "Greys")

high.1n.1r<-unlist(sapply(domain.list, function(x) x$high_1n_r1[x$valid]))
low.1n.1r<-unlist(sapply(domain.list, function(x) x$low_1n_r1[x$valid]))


gc.domains<-unlist(sapply(domain.list, function(x) x$GC_content[x$valid]))

fc.domain.1n.1r<-high.1n.1r/low.1n.1r

#change directory
dir.creat("plots_GCcorrelation")
setwd("plots_GCcorrelation")


plot.fct<-function(y_axes, x_axes, x.name="GC content", y.name="FC MNase high/low", name ){

        pdf(paste( name,".pdf", sep="_"))
          heatscatter(x=x_axes, y=y_axes,colpal=greys[2:9], cor=T,
                cexplot=1, ylab=y.name, xlab=x.name)
        dev.off()

        #regression

        fit<-lm(y_axes~x_axes)

        sink(paste(name,"fit_regression.txt", sep="_"))
          print(summary(fit))
        sink()

        ### filter cooks distance
        get_CD<-function(fit.fkt=fit, cond, filter=0.001){
                select= "cooks_distance"
                cd<-cooks.distance(fit.fkt)
                pdf(paste(cond, select, "distribution.pdf", sep="_"))
                  plot(density(cd), xlim=c(0,0.01), main="cooks_distance")
                  abline(v=filter)
                dev.off()

                idx<-cd< filter
                return(idx)
        }

        sel.cd<-get_CD(cond=name)

        pdf(paste( name,"cd.pdf", sep="_"))
          heatscatter(x=x_axes[sel.cd], y=y_axes[sel.cd],colpal=greys[2:9], cor=T,
                    cexplot=1, ylab=y.name, xlab=x.name)
        dev.off()


        fit.cd<-lm(y_axes[sel.cd]~x_axes[sel.cd])

        sink(paste(name,"fit_regression_cd.txt", sep="_"))
                print(summary(fit.cd))
        sink()

}

plot.fct(y_axes=fc.domain.1n.1r,x_axes=gc.domains,x.name="GC content",
         y.name="FC MNase high/low",name="FC_vs_GC_r1_n1")



#log_FC
plot.fct(y_axes=log2(fc.domain.1n.1r),x_axes=gc.domains,x.name="GC content",
         y.name="log-FC MNase high/low",name="logFC_vs_GC_r1_n1")
