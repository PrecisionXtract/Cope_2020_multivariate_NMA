# plotting survival based on model parameters,  parameters are transformed to be additive
# parameters need to be transformed back to natural scale
# RE example outputs

#CALL LIBRARIES
library(readxl)
library(openxlsx)
library(dplyr)
library(rjags)
library(flexsurv)

#clear workspace 
rm(list = ls())


workfolder = "./"
datafolder = paste(workfolder,"Data/", sep="")
bugsdatafolder = paste(workfolder,"Data/bugsdata/", sep="")
outfolder=paste(workfolder,"Output/NMA_results/", sep="")
figfolder=paste(outfolder, 'NMA survival plots/',sep='')
dir.create(figfolder )

figfolder1=paste(outfolder, 'NMA HR plots/',sep='')
dir.create(figfolder1 )


#LOAD DATA
modellist=c('weibull', 'lognormal', 'llogis', 'gompertz')

outcome_var ='OS'

if (outcome_var=='OS') {outcome_var_name='Overall survival'}

for (i in 1:length(modellist)){
  
  # RE ...........................................................................................
  indataset = paste(outfolder, modellist[i]," NMA RE.xlsx",sep="")
  dataset2=as.data.frame(read_excel(indataset, sheet = 'coda'))
  
  trt=as.data.frame(read_excel(indataset, sheet = 'trt'))
  names(trt)[1]='trt'
  trt$trt=gsub('non ', "non-", trt$trt)
  trt$trt=gsub('Non ', "Non-", trt$trt)
  trt$trt=gsub(' ', " + ", trt$trt)
  
  
  
  
  time=seq(1,96,1)
  maxt=max(time)
  
  comparator_params=list()
  comparator_surv=list()
  comparator_surv_lower=list()
  comparator_surv_upper=list()
  comparator_haz=list()
  hazv=list()
  hr=hr_med=hr_lower=hr_upper=list()
  max_hr_or=1
  
  
  alphas0=dataset2[,grep('alpha',colnames(dataset2))]
  # start looping through treatments
  for (k in 1:nrow(trt)){
    print(paste('start ', trt$trt[k], k))
    alphas=NULL
    alphas=alphas0[,grep(paste0(',',k), colnames(alphas0)),]
    comparator_params[[k]]=apply(alphas,2,median)
    surv1=matrix(NA, nrow=nrow(alphas),ncol=maxt)
    hazv[[k]]=matrix(NA, nrow=nrow(alphas),ncol=length(time))
    
    #### weibull
    if (modellist[i]=='weibull') {
      # point estimates
      a=comparator_params[[k]][2]+1
      b= exp(- log(exp(comparator_params[[k]][1])/a)/a)  
      comparator_surv[[k]]=1-pweibull(time, shape = a, scale = b)
      comparator_haz[[k]]=hweibull(time, shape = a, scale = b)
      
      # CrI from CODA
      acoda=alphas[,2]+1
      bcoda= exp(- log(exp(alphas[,1])/a)/a) 
      codause=data.frame(a=acoda,b=bcoda)
      surv1=t(1-apply(codause, 1, function(yourData) pweibull(q=time, shape = yourData["a"], scale = yourData["b"] )) )
      hazv[[k]]=t(apply(codause, 1, function(yourData) hweibull(time, shape = yourData["a"], scale = yourData["b"] )) )
    } else # end if weibull
      
      #### lognormal
      if (modellist[i]=='lognormal') {
        # point estimates
        a=-comparator_params[[k]][1]/comparator_params[[k]][2]
        b= -1/comparator_params[[k]][2]  
        comparator_surv[[k]]=1-plnorm(time, meanlog = a, sdlog = b)
        comparator_haz[[k]]=hlnorm(time, meanlog = a, sdlog = b)
        
        # CrI from CODA
        acoda=-alphas[,1]/alphas[,2]
        bcoda= -1/alphas[,2]
        codause=data.frame(a=acoda,b=bcoda)
        surv1=t(1-apply(codause, 1, function(yourData) plnorm(q=time, meanlog = yourData["a"], sdlog = yourData["b"] )) )
        hazv[[k]]=t(apply(codause, 1, function(yourData) hlnorm(time, meanlog = yourData["a"], sdlog = yourData["b"] )) )
      }  else # end if lognormal
        
        #### gompertz
        if (modellist[i]=='gompertz') {
          # point estimates
          a=comparator_params[[k]][2]
          b= exp(comparator_params[[k]][1])
          comparator_surv[[k]]=1-pgompertz(time, shape = a, rate = b)
          comparator_haz[[k]]=hgompertz(time, shape = a, rate = b)
          
          # CrI from CODA
          acoda=alphas[,2]
          bcoda= exp(alphas[,1])
          codause=data.frame(a=acoda,b=bcoda)
          surv1=t(1-apply(codause, 1, function(yourData) pgompertz(q=time, shape = yourData["a"], rate = yourData["b"] )) )
          hazv[[k]]=t(apply(codause, 1, function(yourData) hgompertz(time, shape = yourData["a"], rate = yourData["b"] )) )
        } else # end if gompertz
          
          #### log logistic
          if (modellist[i]=='llogis') {
            # point estimates
            a=exp(comparator_params[[k]][2])
            b= 1/exp(comparator_params[[k]][1])
            comparator_surv[[k]]=1-pllogis(time, shape = a, scale = b)
            comparator_haz[[k]]=hllogis(time, shape = a, scale = b)
            
            # CrI from CODA
            acoda=exp(alphas[,2])
            bcoda= 1/exp(alphas[,1])
            codause=data.frame(a=acoda,b=bcoda)
            surv1=t(1-apply(codause, 1, function(yourData) pllogis(q=time, shape = yourData["a"], scale = yourData["b"] )) )
            hazv[[k]]=t(apply(codause, 1, function(yourData) hllogis(time, shape = yourData["a"], scale = yourData["b"] )) )
          }  # end if llogis
            
            
    comparator_surv_lower[[k]]=apply(surv1,2,quantile,probs = c(0.025))
    comparator_surv_upper[[k]]=apply(surv1,2,quantile,probs = c(0.975))
    
    ### hr of treatment vs reference
    hr[[k]]=hazv[[k]]/hazv[[1]]
    hr_med[[k]]=apply(hr[[k]],2,quantile,probs = c(0.5))
    hr_lower[[k]]=apply(hr[[k]],2,quantile,probs = c(0.025))
    hr_upper[[k]]=apply(hr[[k]],2,quantile,probs = c(0.975))
    max_hr_or=max(hr_upper[[k]],max_hr_or)
    
  } # end k loop of treatments
  
  
  
  # max time for plotting, intervals of 3, 6 or 12 months
  if (maxt<=24) {mmgap=3} else {
    if (maxt<=48) {mmgap=6} else{
      mmgap=12
    }
  }
  maxt=ceiling(maxt/mmgap)*mmgap
  
  colour=c('black','green','blue','red')
  
  pngout=paste(figfolder, modellist[i]," NMA survival RE95.png", sep="")
  png(filename=pngout,width=8, height=5, units="in", res=300)
  plot(0,0,type="n",xlab="Time (Months)",ylab=outcome_var_name,
       ylim=c(0,1),xlim=c(0,maxt), las=1,lwd=1.5,yaxt='n',xaxt='n',yaxs="i",xaxs="i", bty='l')
  axis(side=2, las=2)
  axis(1, at=seq(0,maxt,mmgap),labels=seq(0,maxt,mmgap), las=0)
  
  for (k in 1:nrow(trt)){
    lines(c(0,time),c(1,comparator_surv[[k]]),col=colour[k],lwd=3)
    lines(c(0,time),c(1,comparator_surv_lower[[k]]),col=colour[k],lwd=1.5,lty=3)
    lines(c(0,time),c(1,comparator_surv_upper[[k]]),col=colour[k],lwd=1.5,lty=3)
    
  }
  
  legend("topright", legend = c(trt$trt), lty = 1,lwd=3, col = colour,bty='n')
  par(xpd=TRUE)
  text(x=0,y=-0.3*1, labels="The dashed lines represent the 95% credible intervals", cex=0.75, adj=0)
  dev.off()
  
  
  
  ymax=min(5,floor(max_hr_or)+1)
  
  pngout=paste(figfolder1, modellist[i]," NMA hr RE95.png", sep="")
  png(filename=pngout,width=8, height=5, units="in", res=300)
  plot(0,0,type="n",xlab="Time (Months)",ylab=paste('Hazard ratio vs DTIC'),
       ylim=c(0,ymax),xlim=c(0,maxt), las=1,lwd=1.5,yaxt='n',xaxt='n',yaxs="i",xaxs="i", bty='l')
  axis(side=2, las=2)
  axis(1, at=seq(0,maxt,mmgap),labels=seq(0,maxt,mmgap), las=0)
  
  for (k in 2:nrow(trt)){
    lines(time,hr_med[[k]],col=colour[k],lwd=3)
    lines(time,hr_lower[[k]],col=colour[k],lwd=1.5,lty=3)
    lines(time,hr_upper[[k]],col=colour[k],lwd=1.5,lty=3)
  }
  
  legend("topright", legend = paste(c(trt$trt)[2:dim(trt)[1]], 'vs', 'DTIC'), lty = 1,lwd=3, col = colour[2:dim(trt)[1]],bty='n')
  par(xpd=TRUE)
  text(x=0,y=-0.3*ymax, labels="The dashed lines represent the 95% credible intervals", cex=0.75, adj=0)
  dev.off()
  
} #### end i loop of distributions

