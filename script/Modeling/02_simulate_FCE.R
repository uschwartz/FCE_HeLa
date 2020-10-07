args<-commandArgs(TRUE)

project.path<-args[1]
working.path<-paste0(project.path,"/modeling/02_digest")

dir.create(working.path, recursive = T)

## working path
setwd(working.path)
dir.create("plots_domain")

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

plotNucDist<-function(ir.data, cond, temp){
        mono<-ir.data[width(ir.data)==181,]
        mono.trim<-trimRange(irange=mono, trim=30)
        all<-mono.trim       
        cov<-as.vector(coverage(all)) 
        pdf(paste(file="plots_domain/",temp,"_",cond,".pdf",sep=""))
        plot(cov,type="l", lty=1, ylim=c(10000,30000),ylab="coverage", main=cond, xlab="Position" )
        dev.off()
}



getMono<-function(cut.freq.in, template.in){
        x<-(-1)
        while(x<1000000){    
                ir.data<-Model.MNase(cut.freq=cut.freq.in, runs=10000, template=template.in)
                if(x==-1){
                        mono<-ir.data[width(ir.data)==181,]
                } else{
                       mono<-c(mono,ir.data[width(ir.data)==181,]) 
                }
                x<-length(mono)
                print(x)
        }
        return(mono[1:1000000])       
}

################################################################################
############### open chromatin configuration ###################################
#################################################################################

numb.nucs<-50
linker.length<-30

cut.pos<- seq(0,((150+linker.length)*(numb.nucs)-90),90)
shift.nuc<-seq(2,100,4)
cut.pos[shift.nuc]<-cut.pos[shift.nuc]-5

prob.cut<-rep(c(10,1), numb.nucs)
standard.stretch<-data.frame(cut_site=cut.pos,cut_prob=prob.cut)

## high-MNase simulation
mono.norm.sU70.std<-getMono(cut.freq.in=70,template.in=standard.stretch)
save(mono.norm.sU70.std, file="mono.norm.sU70.std.rda")

plotNucDist(ir.data=mono.norm.sU70.std, cond="70sU", temp="standard")

## low-MNase simulation
mono.norm.sU20.std<-getMono(cut.freq.in=20,template.in=standard.stretch)
save(mono.norm.sU20.std, file="mono.norm.sU20.std.rda")

plotNucDist(ir.data=mono.norm.sU20.std, cond="20sU", temp="standard")


################################################################################
############### condensed chromatin domain in the center #######################
#################################################################################

prob.cut.normal<-rep(c(10,1), numb.nucs/2)
prob.cut.cond<-rep(c(1,8), numb.nucs/2)

cond.stretch<-data.frame(cut_site=cut.pos,cut_prob=c(prob.cut.normal[1:25],prob.cut.cond,
                prob.cut.normal[26:50]))
#low-MNase
mono.norm.sU20.cond<-getMono(cut.freq.in=20,template.in=cond.stretch)
save(mono.norm.sU20.cond, file="mono.norm.sU20.cond.rda")
plotNucDist(ir.data=mono.norm.sU20.cond, cond="20sU", temp="condensed_8")

#high-MNase
mono.norm.sU70.cond<-getMono(cut.freq.in=70,template.in=cond.stretch)
save(mono.norm.sU70.cond, file="mono.norm.sU70.cond.rda")
plotNucDist(ir.data=mono.norm.sU70.cond, cond="70sU", temp="condensed_8")



################################################################################
############### AT rich domain in the center #######################
#################################################################################

# change accessibility of nucleosome core

#standard stretch
numb.nucs<-50
linker.length<-30

for (k in c(1.5,2.5,3)){
    prob.cut.normal<-rep(c(10,1)*2, numb.nucs/2)
    prob.cut.inst<-rep(c(k,10)*2, numb.nucs/2)
    
    instabil.stretch<-data.frame(cut_site=cut.pos,
                cut_prob=c(prob.cut.normal[1:25],prob.cut.inst,
                                    prob.cut.normal[26:50]))
    
    mono.norm.sU20.inst<-getMono(cut.freq.in=20,template.in=instabil.stretch)
    mono.norm.sU70.inst<-getMono(cut.freq.in=70,template.in=instabil.stretch)
    
    save(mono.norm.sU20.inst, file="mono.norm.sU20.inst.rda")
    save(mono.norm.sU70.inst, file="mono.norm.sU70.inst.rda")
    
    
    plotNucDist(ir.data=mono.norm.sU20.inst, cond="20sU",
                temp=paste0("instabil_stretch_",k))
    plotNucDist(ir.data=mono.norm.sU70.inst, cond="70sU",
                temp=paste0("instabil_stretch_",k))
    
}




