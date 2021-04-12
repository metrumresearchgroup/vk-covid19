# Script name: COVID_viral_load_dose_simulation.R
# Analyst name: Emmanuel Chigutsa
# Date created: April 2020
# Software version: R version 3.5.0 for Windows

library(RxODE)
library(scales)
library(ggplot2)
require(gridExtra)
library(truncnorm)
library(dplyr)
library(gridExtra)
rm(list=ls())

####### define model ######
ode<-"
  KD=ec50;
  conc=centr/vc;  # drug concentration in nM
  lungconc=delay*scalar;  # lung concentration is x% that of serum for mAbs  
 fb = lungconc/(lungconc+KD) ; # fraction virus bound, KD in nM

 virusfree=VL*(1-fb);  # virus not bound to LY, available for infecting target cells

   d/dt(centr) = -K10*centr - K12*centr + K21*peri;  # mAb central compt
   d/dt(delay) = keo*(conc-delay);
   d/dt(peri) =  K12*centr - K21*peri;  # mAb peripheral cmpt
   d/dt(dft) = - beta*dft*virusfree;  # beta is rate constant of infection
   d/dt(VL)=gamma*dft*virusfree - delta*VL;  # total virus
   d/dt(vauc) = VL;  # viral load auc
"
##### end of model ####

mod4<-RxODE(model = ode, modName = "mod4")

# list model parameters
npat=12000  #12000 patients

# need function to chop off 5th and 95th tails #
chop<-function(x){
   x<-x[x>(quantile(x,0.05))&x<(quantile(x,0.95))] 
}

##### viral kinetics
# VL0 = 21.8*exp(rnorm(n=npat, mean=0, sd=3.65))
VL0 = 21.8*exp(rnorm(n=npat, mean=0, sd=2))
VL0 = chop(VL0)
beta=(6.77e-5)*exp(rnorm(n=npat, mean=0, sd=1.75))/5  # infection rate constant ((copies/mL)-1)/day)
beta=chop(beta)
gamma=3.8*exp(rnorm(n=npat, mean=0, sd=0.0787))  # maximum rate constant for virus infection
gamma=chop(gamma)
delta=0.7*exp(rnorm(n=npat, mean=0, sd=0.0464))  # (/day)  # death rate of infected cells
delta=chop(delta)
#######

ec50=1.1  # 

## PK ### sourced from canakinumab paper, IgG1 antibody. Chakraborty et al. 2012.
cl = 0.137*exp(rnorm(n=npat, mean=0, sd=0.3))  # (L/d/70 kg)
cl=chop(cl)
vc = 3.3*exp(rnorm(n=npat, mean=0, sd=0.3))  # L/70 kg
vc=chop(vc)
vp = 2.71*exp(rnorm(n=npat, mean=0, sd=0.3)) # L/70 kg
vp=chop(vp)
q = 0.429 # L/day/70 kg
K10 = cl/vc
K12 = q/vc
K21 = q/vp
scalar = 0.15
keo=41 # Assume 2 hour for SS lung distribution from serum. 2 h is thalf*5. keo=0.693/thalf, convert to days
dstart=runif(npat,0.5,3)
dstart=chop(dstart)


params.all<-cbind.data.frame(VL0=VL0,beta=beta,gamma=gamma,delta=delta,
                         vc=vc, cl=cl, vp=vp, q=q, K10=K10, K12=K12, K21=K21,ec50=ec50,scalar=scalar,keo=keo,dstart=dstart)
params.all$ID<-1:length(params.all$dstart)
doses=c(0,700/70,1400/70)  # mg/kg
wt=70   # typical weight
dft0=1
params.all<-params.all[1:10000,]

### day 0 ################
res <- NULL #Create an empty matrix for storing results
df.full<-NULL

for(jj in unique(doses)){
   amt=jj
df_jj <- do.call(rbind, lapply(1:nrow(params.all), FUN = function(uu) {
   ev <- eventTable()  
  # Specify sampling
   ev$add.sampling(c(seq(0,7,0.5),seq(8,28,1)))
   params <- params.all[uu,]
    inits<-c(centr=0,delay=0,peri=0,dft=dft0,infected=0,VL=params.all[uu,"VL0"],AUCV=0)  # define initial conditions
    ev$add.dosing(dose=1000000*amt*wt/150000,nbr.doses=1,dosing.to=1,start.time=params.all[uu,"dstart"])  # dose in mg, times 1 million = ng, divide by 150 kDa to get nanomoles
    x <- mod4$run(params, ev, inits = inits)
    x<-cbind(x,"ID"=params.all[uu,"ID"],"ARM"=jj,"START"=params.all[uu,"dstart"])    # append results
    res <- rbind.data.frame(res,x)  # append results into dataframe
})
)
  df.full<-rbind(df.full,df_jj)
}
res<-df.full

data<-data.frame(res)
data$VL <- data$VL*(1 + rnorm(length(data$VL),mean=0,sd=0.2))  # add 20% residual error
names(data)

data.out<-cbind.data.frame(ID=data$ID,ARM=ifelse(data$ARM==0,"Control",
 							ifelse(data$ARM==20,"LY dose 1","LY dose 2")),
                            TIME=data$time,VL=data$VL,VAUC=data$vauc,RAND_START=data$START)

# write.csv(data.out,"mild_moderate_patients_standard_variability_2.csv",quote=F, na=".",row.names=F)


data.plot<-data%>%
  group_by(time,ARM)%>%
  summarise(
    lowpi=quantile(VL,0.025),
    med=quantile(VL,0.5),
    highpi=quantile(VL,0.975)
  )

d.p<-cbind.data.frame(xx=c(0.5,3),ymin=10e-6,ymax=10e6)

pdf("virtual_patients_plots_700_1400.pdf",width=7,height=6)
ggplot(data.plot,aes(x=time, y=med,fill=as.factor(ARM))) +
 geom_ribbon(data=d.p,aes(x=xx,ymin=ymin,ymax=ymax),alpha=0.3,inherit.aes = FALSE,fill='grey4')+
  geom_ribbon(aes(x=time,ymin=lowpi,ymax=highpi),alpha=0.3)+
geom_line(cex=1.2,aes(color=factor(ARM)))+
geom_hline(yintercept = 1,col="black",cex=1.2,lty="dashed")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position=c(0.8,0.8),legend.text=element_text(size=12),legend.title=element_text(size=12))+
  scale_x_continuous(breaks=seq(0,28,7),limits=c(0,28),labels=seq(0,28,7))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)),limits=c(0.00001,10000000)) +
 # ggtitle("Viral load over time - 95% prediction interval \nMild-moderate patients (standard variability)") +
  xlab("Time since onset of symptoms (days)") + 
  ylab("Viral load (copies/mL)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title=element_text(size=17))+
scale_fill_discrete(name="Study arms",
                        labels=c("Control", "700 mg","1400 mg"))+
scale_color_discrete(name="Study arms",
                        labels=c("Control", "700 mg","1400 mg"))
dev.off()


