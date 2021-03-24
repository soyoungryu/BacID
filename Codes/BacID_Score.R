directory = '[SPECIFY DIRECTORY]'
setwd(directory)
source('bacID_functions.R')
load('Data.RData')

#PARAMETERS
no.cores=500 # No. of available cores for this analysis

#Reference spectra mass recalibration. 
#Mass recalibration values were obtained from a train set
mz.shift = c(0, -3.251  , -2.379   , -3.915   , 0.328   , 0.124, 0, -2.147 , rep(0, length(ref)-8)) 
ref = mass.recal(ref, mz.shift) #ref is already normalized. No need to do it again. 
#Normalize intensity values of observed spectra
Mix$normalized.peak.list <- list()
for (i in 1:length(Mix$peak.list)){
Mix$normalized.peak.list[[i]] = normalize.intensity(Mix$peak.list[[i]], log.transform=TRUE)
}

#Fit Model & Calculate BacID Prob
library(doParallel)
cl <- makeCluster(no.cores)
registerDoParallel(cl)
results <- foreach(i1=1:length(Mix$normalized.peak.list)) %dopar% {
  poly <- Mix$normalized.peak.list[[i1]][,1:2]; colnames(poly) <- c('mz', 'int')
  true.id <- colnames(Mix$microbe.list)[Mix$microbe.list[i1,]]
  conf_score(poly, ref, true.id)
}
stopCluster(cl)
save(results, file="BacID_Score_Results.RData")

results <- results[Mix$dataset.name=="Standard"]
log.probs = id = correct = file.name = NULL
for (i in 1:length(results)){
log.probs <- c(log.probs, results[[i]]$log.probs)
id <-  c(id, results[[i]]$id) #bacterial IDs
correct <- c(correct, results[[i]]$correct.id)
}

pvalue <- p_estimate(log.probs, correct) #estimate p-values using decoy results
fdr <- p.adjust(pvalue, method="fdr") #remove decoy results

sum(fdr[correct>=0]<0.01) #no. of bacteria with fdr<0.01 among non-decoy results (correct=0; incorrect id, correct=1 for correct id, and correct=-1 for decoy id)
sum(correct[fdr<0.01]==0)/sum(fdr[correct>=0]<0.01) #true fdr 





