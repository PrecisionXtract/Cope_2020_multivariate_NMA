## This file uses R relative paths. Open this file through R project file "Cope_2020_multivariate_NMA.Rproj" in the top "Cope_2020_multivariate_NMA" directory
## running multivariate normal NMA in jags
## examples using weibull, log norm, log logis, and gompertz parameters

#CALL LIBRARIES
library(readxl)
library(openxlsx)
library(dplyr)
library(rjags)

#clear workspace 
rm(list = ls())


workfolder = "./"
datafolder = paste(workfolder,"Data/", sep="")
bugsdatafolder = paste(workfolder,"Data/bugsdata/", sep="")
outfolder=paste(workfolder,"Output/NMA_results/", sep="")
jagscodefolder=paste(workfolder,"Jagscode/", sep="")



#LOAD DATA
modellist=c('weibull', 'lognormal', 'llogis', 'gompertz')

for (i in 1:length(modellist)){
  model1=modellist[i]
  
  dataname=paste(bugsdatafolder, 'parametric survival bugs data.xlsx', sep="")
  
  dataset2=as.data.frame(read_excel(dataname, sheet=paste('data2',model1)))
  
  dataset3=as.data.frame(read_excel(dataname, sheet="data3"))
  
  trt=as.data.frame(read_excel(dataname, sheet="treatments"))
  
  
  ## prepping data to run in bugs  
  N1=dim(dataset2)[1]              #no of datapoints
  N2=dim(dataset3)[1]             # no of studies x no of outcomes (10x2=20) 
  ns=length(unique(dataset3$studyid))              # no of studies
  no=num_o= dim(dataset2[, grep('y', substr(colnames(dataset2),1,1))])[2]      #no of outcomes
  nt.total =rep(length(unique(dataset2$treatment)),no)                     # no of interventions for each outcome 
  nastudy = c(table(dataset2$studyid)) # no. of arms in each study
  
  studyid=dataset2$studyid
  study=dataset2$study
  y=as.matrix(dataset2[,grep('y\\.',colnames(dataset2))])
  se=as.matrix(dataset2[,grep('se\\.',colnames(dataset2))])
  arm=dataset2$arm
  cov11= dataset2$cov11
  cov22= dataset2$cov22
  cov12= dataset2$cov12
  
  
  studyid1= dataset3$studyid1
  t=cbind(dataset3$t1,dataset3$t2,dataset3$t3)
  s= dataset3$s	
  o= dataset3$o	
  out= dataset3$out	
  na= dataset3$na
  
  
  
  datause = list(    N1      = N1,
                     N2      = N2         ,
                     ns      = ns         ,
                     no      = no         ,
                     nt.total= nt.total   ,
                     nastudy = nastudy    ,
                     studyid = as.numeric(as.factor(studyid))    ,
                     study   = study      ,
                     arm     = arm        ,
                     y       = y          ,
                     cov11   = cov11      ,
                     cov22   = cov22      ,
                     cov12   = cov12      ,
                     studyid1= as.numeric(as.factor(studyid1))  ,
                     s       = s          ,
                     t       = t          ,
                     o       = o          ,
                     out     = out        ,
                     na      = na         
  )
  
  
  
  #NUMBER OF ITERATIONS
  iterate=20000 
  NTHIN=5
  iterate1=iterate*NTHIN
  
  
  # Fixed effect model
  MODELFILE<-paste(jagscodefolder, "FE 2outcomes/Anchana 2014 code FE.bugs", sep="" ) 
  set.seed(1)
  
  jagsFE=NULL
  jagsFE<- suppressWarnings(jags.model(file = MODELFILE, data =datause, n.chains = 2, n.adapt = 80000))
  update(jagsFE,iterate)
  jagsFE.out <- coda.samples(jagsFE, c("d","mu", 'mu_mean','mean.y','alpha'), n.iter = iterate1, thin = NTHIN)
  DICFE = dic.samples(jagsFE, iterate1, thin = NTHIN, type= "pD")
  summaryFE = summary(jagsFE.out)
  rawFE = data.frame(variable=rownames(summaryFE$statistics),(summaryFE$statistics)[,c(1,2)],summaryFE$quantiles, summaryFE$start, summaryFE$end)
  pD = round(sum(DICFE$penalty),2)
  deviance = round(sum(DICFE$deviance),2)
  DIC0 =  pD + deviance
  DIC=data.frame(Dbar=deviance, pD=pD, DIC=DIC0)
  
  # coda
  codadataset3=do.call(rbind.data.frame, mcmc.list(jagsFE.out[[1]],jagsFE.out[[2]]))
  dim(codadataset3)
  print(paste('total coda sample, all chains = ', dim(codadataset3)[1]))
  
  codakeep.FE=codadataset3[,c(grep("d\\[",substr(colnames(codadataset3),1,2) ),grep("alpha\\[",colnames(codadataset3)) ,grep("mu_mean",colnames(codadataset3)) ) ]
  colnames(codakeep.FE)
  
  
  outputfile = paste(outfolder, model1," NMA FE.xlsx",sep="")
  
  outdata=list(rawFE, trt, DIC, dataset2,codakeep.FE)
  openxlsx::write.xlsx(outdata,outputfile,
                       sheetName= c('raw bugs output FE','trt', 'DIC','bugs data', 'coda' ) )
  
  
  
  # Random effects model
  MODELFILE<-paste(jagscodefolder, "RE 2outcomes/Anchana 2014 code RE.bugs", sep="" ) 
  set.seed(1)
  
  jagsRE=NULL
  jagsRE<- suppressWarnings(jags.model(file = MODELFILE, data =datause, n.chains = 2, n.adapt = iterate))
  update(jagsRE,iterate)
  jagsRE.out <- coda.samples(jagsRE, c("d","mu", 'mu_mean','mean.y','alpha'), n.iter = iterate1, thin = NTHIN)
  DICRE = dic.samples(jagsRE, iterate1, thin = NTHIN, type= "pD")
  summaryRE = summary(jagsRE.out)
  rawRE = data.frame(variable=rownames(summaryRE$statistics),(summaryRE$statistics)[,c(1,2)],summaryRE$quantiles, summaryRE$start, summaryRE$end)
  pD = round(sum(DICRE$penalty),2)
  deviance = round(sum(DICRE$deviance),2)
  DIC0 =  pD + deviance
  DIC=data.frame(Dbar=deviance, pD=pD, DIC=DIC0)
  
  
  # coda
  codadataset3=do.call(rbind.data.frame, mcmc.list(jagsRE.out[[1]],jagsRE.out[[2]]))
  dim(codadataset3)
  print(paste('total coda sample, all chains = ', dim(codadataset3)[1]))
  
  codakeep.RE=codadataset3[,c(grep("d\\[",substr(colnames(codadataset3),1,2) ),grep("alpha\\[",colnames(codadataset3)) ,grep("mu_mean",colnames(codadataset3)) ) ]
  colnames(codakeep.RE)
  
  
  outputfile = paste(outfolder, model1," NMA RE.xlsx",sep="")
  
  outdata=list(rawRE, trt, DIC, dataset2,codakeep.RE)
  openxlsx::write.xlsx(outdata,outputfile,
                       sheetName= c('raw bugs output RE','trt', 'DIC','bugs data','coda' ) )
  
} # end looping of models



## Next step is to take alphas from raw bugs output files, transform back to appropriate scale for given suvival function, and create survival curves, hazard curves etc.
