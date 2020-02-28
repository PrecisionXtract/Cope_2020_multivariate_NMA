## This file uses R relative paths. Open this file through R project file "Cope_2020_multivariate_NMA.Rproj" in the top "Cope_2020_multivariate_NMA" directory
## jags needs to be downloaded and installed on your computer http://mcmc-jags.sourceforge.net/

## INSTALL THE FOLLOWING LIBRARIES
if(!require(readxl)) {install.packages("readxl")}
if(!require(openxlsx)) {install.packages("openxlsx")}
if(!require(survival)) {install.packages("survival")}
if(!require(flexsurv)) {install.packages("flexsurv")}
if(!require(dplyr)) {install.packages("dplyr")}
if(!require(rjags)) {install.packages("rjags")}


## check that they are installed
'readxl' %in% rownames(installed.packages()) 
'openxlsx' %in% rownames(installed.packages()) 
'survival' %in% rownames(installed.packages()) 
'flexsurv' %in% rownames(installed.packages()) 
'dplyr' %in% rownames(installed.packages()) 
'rjags' %in% rownames(installed.packages()) 


