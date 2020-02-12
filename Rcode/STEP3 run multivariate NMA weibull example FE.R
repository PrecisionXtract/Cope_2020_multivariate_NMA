## running multivariate NMA
## example using weibull parameters


#CALL LIBRARIES
library(readxl)
library(openxlsx)
library(dplyr)
library(R2WinBUGS)
library("coda")

#clear workspace 
rm(list = ls())


workfolder = "./"
datafolder = paste(workfolder,"Data/", sep="")
bugsdatafolder = paste(workfolder,"Data/bugsdata/", sep="")

## winbugs does not function with relative pathways
wd=getwd()
outfolder=paste(wd,"/Output/NMA_results/", sep="")
bugscodefolder=paste(wd,"/Bugscode/", sep="")
bugs_directory = paste(wd,"/WinBUGS14/", sep="")


dir.create(outfolder )



#LOAD DATA
model1=c('weibull')

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
R2WinBUGS::bugs.data(datause)



#NUMBER OF ITERATIONS
NTHIN<-2
NBURNIN<-20000*NTHIN
NITER<-2*NBURNIN


# Fixed effect model
MODELFILE<-paste(bugscodefolder, "FE 2outcomes/Anchana 2014 code FE.bugs", sep="" ) 
output_location=outfolder

NMA.sim=NULL
NMA.sim<- R2WinBUGS::bugs(data=datause, inits=NULL, model.file = MODELFILE,
                           parameters=c("d","mu", 'mu_mean','mean.y','alpha'),
                           n.chains = 2, n.iter = NITER, n.thin= NTHIN, DIC=TRUE, debug=FALSE, 
                           bugs.seed=1, codaPkg = FALSE, working.directory=output_location, 
                          save.history=FALSE,bugs.directory = bugs_directory)

rawFE = NMA.sim$summary
pD = round(NMA.sim$pD,2)
NMA.sim$DIC
DIC=data.frame(Dbar=NMA.sim$DIC-pD, pD=pD, DIC=NMA.sim$DIC)
unlink(paste0(output_location,"coda1.txt"))
unlink(paste0(output_location,"coda2.txt"))
unlink(paste0(output_location,"codaIndex.txt"))
unlink(paste0(output_location,"data.txt"))
unlink(paste0(output_location,"log.txt"))
unlink(paste0(output_location,"log.odc"))
unlink(paste0(output_location,"script.txt"))

outputfile = paste(outfolder, model1," NMA FE.xlsx",sep="")

outdata=list(rawFE, trt, DIC, dataset2)
openxlsx::write.xlsx(outdata,outputfile,
                     sheetName= c('raw bugs output FE','trt', 'DIC','bugs data' ),rowNames=c(TRUE,FALSE,FALSE,FALSE) )



# Random effects model
MODELFILE<-paste(bugscodefolder, "RE 2outcomes/Anchana 2014 code RE.bugs", sep="" ) 
output_location=outfolder

NMA.sim=NULL
NMA.sim<- R2WinBUGS::bugs(data=datause, inits=NULL, model.file = MODELFILE,
                          parameters=c("d","mu", 'mu_mean','mean.y','alpha'),
                          n.chains = 2, n.iter = NITER, n.thin= NTHIN, DIC=TRUE, debug=FALSE, 
                          bugs.seed=1, codaPkg = FALSE, working.directory=output_location, 
                          save.history=FALSE,bugs.directory = bugs_directory)

rawRE = NMA.sim$summary
pD = round(NMA.sim$pD,2)
NMA.sim$DIC
DIC=data.frame(Dbar=NMA.sim$DIC-pD, pD=pD, DIC=NMA.sim$DIC)
unlink(paste0(output_location,"coda1.txt"))
unlink(paste0(output_location,"coda2.txt"))
unlink(paste0(output_location,"codaIndex.txt"))
unlink(paste0(output_location,"data.txt"))
unlink(paste0(output_location,"log.txt"))
unlink(paste0(output_location,"log.odc"))
unlink(paste0(output_location,"script.txt"))

outputfile = paste(outfolder, model1," NMA RE.xlsx",sep="")

outdata=list(rawRE, trt, DIC, dataset2)
openxlsx::write.xlsx(outdata,outputfile,
                     sheetName= c('raw bugs output RE','trt', 'DIC','bugs data' ),rowNames=c(TRUE,FALSE,FALSE,FALSE) )


## Next step is to take alphas from raw bugs output files, transform back to appropriate scale for given suvival function, and crate survival curves
