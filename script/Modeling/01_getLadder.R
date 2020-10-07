args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-paste0(project.path,"/modeling/01_ladder")

dir.create(working.path, recursive = T)

## working path
setwd(working.path)

library(IRanges)


Model.MNase<-function(cut.freq, runs, template){
        #create cut.prob
        cut.prob<-c()
        for( i in 1:nrow(template)){
                cut.prob<-c(rep(i, template$cut_prob[i]),cut.prob)
        }
        frags<-1
        # occupy possible restriction sites
        for(j in 1:runs){
                drawn<-c()
                new.prob<-cut.prob
                for(i in 1:cut.freq){
                        draw<-sample(new.prob, 1)
                        drawn<-c(drawn,draw)
                        new.prob<-new.prob[new.prob!=draw]
                }
                #fragments
                cuts<-template$cut_site[drawn]       
                cuts.ord<-cuts[order(cuts)]
                ir.frags<-IRanges(start=cuts.ord[1:(cut.freq-1)], end=cuts.ord[2:cut.freq])
                
                if(length(frags)==1)  frags<-ir.frags else frags<-c(frags,ir.frags)
                print(j)
        }
        return(frags)
}


trimRange<-function(irange,trim){
        start(irange)<- start(irange)+trim
        end(irange)<- end(irange)-trim
        return(irange)
}

#standard stretch
numb.nucs<-50
linker.length<-30

cut.pos<- seq(0,((150+linker.length)*(numb.nucs)-90),90)
shift.nuc<-seq(2,100,4)
cut.pos[shift.nuc]<-cut.pos[shift.nuc]-5


prob.cut<-rep(c(10,1), numb.nucs)
standard.stretch<-data.frame(cut_site=cut.pos,cut_prob=prob.cut)


ir.sU20.std<-Model.MNase(cut.freq=20,runs=10000,template=standard.stretch)
save(ir.sU20.std, file="ir.sU20.std.rda")

ir.sU70.std<-Model.MNase(cut.freq=70,runs=10000,template=standard.stretch)
save(ir.sU70.std, file="ir.sU70.std.rda")




pdf("raw_20sU_ladder.pdf")
    hist(log10(width(ir.sU20.std)), main="20 sU - MNase", col="green",
         xlim=c(1.7,log10(10000)), xaxt="n", xlab="fragment length")
    axis(1,at=log10(c(180, 380, 530, 730, 1030)), 
         labels=as.character(c(150, 350, 500, 700, 1000)))
dev.off()


pdf("raw_70sU_ladder.pdf")
    hist(log10(width(ir.sU70.std)), main="70 sU - MNase", col="darkred"  ,
       xlim=c(1.7,log10(10000)), xaxt="n", xlab="fragment length")
    axis(1,at=log10(c(180, 380, 530, 730, 1030)), 
       labels=as.character(c(150, 350, 500, 700, 1000)))
dev.off()




# Monos Convert to gel intensity / log scale
getGelIntensity<-function(iranges_frags){
        frags.tab<-table(width(iranges_frags))
        intensities<-c()
        for(i in 1:length(frags.tab)){
                frag.length<-as.numeric(names(frags.tab[i]))
                normFrags<-rep(frag.length,frag.length*frags.tab[i])
                intensities<-c(intensities, normFrags)
        }
        return(intensities) 
}


hist.U20<-getGelIntensity(ir.sU20.std)
hist.U70<-getGelIntensity(ir.sU70.std)



pdf("length_Norm_20sU_ladder.pdf")
    hist(log10(hist.U20), main="20 sU - MNase", col="green",breaks=20  ,
     xlim=c(1.7,log10(10000)), 
     xaxt="n", xlab="fragment length", ylim=c(0,6e7))
    axis(1,at=log10(c(180, 380, 530, 730, 1030)), 
     labels=as.character(c(150, 350, 500, 700, 1000)))
dev.off()



pdf("length_Norm_70sU_ladder.pdf")
    hist(log10(hist.U70), main="70 sU - MNase", col="darkred",breaks=10 ,
       xlim=c(1.7,log10(10000)), xaxt="n", xlab="fragment length",
       ylim=c(0,6e7))
    axis(1,at=log10(c(180, 380, 530, 730, 1030)), 
       labels=as.character(c(150, 350, 500, 700, 1000)))
dev.off()



