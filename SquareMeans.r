require(tidyr)
require(arm)

n.ind<-800
n.sim=100
Vs.error<-c(7, 1.28571) 
Vs.trait<-c(3)

mean.trait<-0


l<-list()
l1<-list()
l2<-list()  
lmm<-list()  

varT21<-data.frame(sim=1:n.sim*2, reap=NA, reap2=NA,   var=NA, type="m2")
varT22<-data.frame(sim=1:n.sim*2, reap=NA,reap2=NA,    var=NA, type="2m")
varT2<-data.frame(sim=1:n.sim*2, reap=NA, reap2=NA,   var=NA, type="real")
varTmm<-data.frame(sim=1:n.sim*2, reap=NA, reap2=NA,   var=NA, type="mm")

for(j in 1:2){
V.trait<-Vs.trait
V.error<-Vs.error[j]
for(i in 1:n.sim){
  d<-data.frame(ID=1:n.ind, t=rnorm(n.ind, mean.trait, sqrt(V.trait)))
  d$t2<-d$t^2
  d$z<-scale(d$t)
  d$z2<-scale(d$t2)
  
  d$obs.t1<-rnorm(n.ind, d$t, sqrt(V.error))
  d$obs.t2<-rnorm(n.ind, d$t, sqrt(V.error))
  d$obs.t3<-rnorm(n.ind, d$t, sqrt(V.error))
  
  d$obs.t12<-d$obs.t1^2
  d$obs.t22<-d$obs.t2^2
  d$obs.t32<-d$obs.t3^2
  
d$mean.t<-apply(d[,c("obs.t1", "obs.t2", "obs.t3")],1,mean)
d$mean.t2<-apply(d[,c("obs.t12", "obs.t22", "obs.t32")],1,mean)
d$mean.t22<-d$mean.t^2

d2.b <- gather(d, obs, t2x, obs.t12:obs.t32, factor_key=TRUE)

mod<-lmer(t2x~1 + (1|ID), data=d2.b)
VarCorr(mod)$ID[1,1]

varTmm$var[i]<-VarCorr(mod)$ID[1,1]
varT2$var[i]<-var(d$t2)
varT21$var[i]<-var(d$mean.t2)
varT22$var[i]<-var(d$mean.t22)

varTmm$reap[i]<-V.trait/(V.trait+V.error)
varT2$reap[i]<-V.trait/(V.trait+V.error)
varT21$reap[i]<-V.trait/(V.trait+V.error)
varT22$reap[i]<-V.trait/(V.trait+V.error)

varTmm$reap2[i]<-var(d$t2)/var(d2.b$t2x)
varT2$reap2[i]<-var(d$t2)/var(d2.b$t2x)
varT21$reap2[i]<-var(d$t2)/var(d2.b$t2x)
varT22$reap2[i]<-var(d$t2)/var(d2.b$t2x)

  
}
l[[j]]<-varT2
l1[[j]]<-varT21
l2[[j]]<-varT22

lmm[[j]]<-varTmm
}


              VT2I<-2*Vs.trait^2 + 4*Vs.trait*mean.trait^2
  VT2R<-2*(Vs.trait +Vs.error)^2 + 4*(Vs.trait+Vs.error)*mean.trait^2
VT2I/(VT2R)



x<-do.call(rbind.data.frame, l)
x1<-do.call(rbind.data.frame, l1)
x2<-do.call(rbind.data.frame, l2)
xmm<-do.call(rbind.data.frame, lmm)

d<-rbind(x1,x2,xmm)
d$reap<-round(d$reap,2)

p<- ggplot(d, aes(x=as.factor(reap), y=var, fill=type)) + 
  geom_boxplot() 
p +labs(title="", x="Repeatability", y = "Among individual variance", fill="Repeatability")+
  theme_classic() +geom_hline(yintercept=VT2I, linetype="dashed") + scale_fill_manual(values=c("white", "darkgrey", "black"), 
                                                                                      name="Models",
                                                                                      breaks=c("m2", "2m", "mm"),
                                                                                      labels=c("mean of squared values", "square of mean value", "mixed model estimate" )) 


plot(log(xmm$var)~log(x$var), col=as.numeric(as.factor(c(x$reap,xmm$reap))))


(-1.5/sqrt(0.3))/(2*(-0.3/(0.3)))
(-1.5)/(2*(-0.3))
        
