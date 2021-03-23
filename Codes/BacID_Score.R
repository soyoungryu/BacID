directory = '[SPECIFY DIRECTORY]'

setwd(directory)
source('bacID_functions.R')
load('Data.RData')

#PARAMETERS
no.cores=500 # No. of available cores for this analysis

#Reference spectra mass recalibration. 
#Mass recalibration values were obtained from a train set
mz.shift = c(0, -3.251  , -2.379   , -3.915   , 0.328   , 0.124, 0, -2.147 , rep(0, length(ref)-8)) 
ref = mass.recal(ref, mz.shift)
for (i in 1:length(ref)){ref[[i]] = normalize.intensity(ref[[i]], log.transform=TRUE)}
#Normalize intensity values of observed spectra
Mix$normalized.peak.list <- list()
for (i in 1:length(Mix$peak.list)){
Mix$normalized.peak.list[[i]] = normalize.intensity(Mix$peak.list[[i]], log.transform=TRUE)
}

#Fit Model & Compute BacID Scores
library(doParallel)
cl <- makeCluster(no.cores)
registerDoParallel(cl)
results <- foreach(i1=1:length(Mix$normalized.peak.list)) %dopar% {
  mass.tolerance=2000
  poly <- Mix$normalized.peak.list[[i1]][,1:2]; colnames(poly) <- c('mz', 'int')
  true.id <- colnames(Mix$microbe.list)[Mix$microbe.list[i1,]]
  conf_score(poly, ref, true.id, mass.tolerance)
}
stopCluster(cl)

save(file='BacID_Score_Results.RData')

results <- results[Mix$dataset.name=="Standard"]
scores = id = correct = file.name = NULL
for (i in 1:length(results)){
scores <- c(scores, results[[i]]$scores)
id <-  c(id, results[[i]]$ids) #bacterial IDs
correct <- c(correct, results[[i]]$correct.id)
}

pvalue <- p_estimate(scores, correct) #estimate p-values using decoy results
fdr <- p.adjust(pvalue, method="fdr") 

#No of IDs at estimated FDRs < 1% among non-decoy identificaitons (correct=0 or 1)
sum(fdr[correct>=0]<=0.01) # 
#True FDR at estimated FDRs of 1%
sum(correct[fdr<=0.01]==0)/sum(fdr[correct>=0]<=0.01) #






