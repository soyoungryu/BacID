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
Mix.test$normalized.peak.list <- list()
for (i in 1:length(Mix.test$peak.list)){
Mix.test$normalized.peak.list[[i]] = normalize.intensity(Mix.test$peak.list[[i]], log.transform=TRUE)
}

#Fit Model & Calculate BacID Prob
library(doParallel)
cl <- makeCluster(no.cores)
registerDoParallel(cl)
results <- foreach(i1=1:length(Mix.test$normalized.peak.list)) %dopar% {
  mass.tolerance=2000
  poly <- Mix.test$normalized.peak.list[[i1]][,1:2]; colnames(poly) <- c('mz', 'int')
  true.id <- colnames(Mix.test$microbe.list)[Mix.test$microbe.list[i1,]]
  matching_bi(poly, ref, true.id, mass.tolerance)
}
stopCluster(cl)

save(file='BacID_Prob_Results.RData')

results <- results[Mix.test$dataset.name=="Standard"]
log.probs = id = correct = file.name = NULL
for (i in 1:length(results)){
log.probs <- c(log.probs, results[[i]]$log.probs)
id <-  c(id, results[[i]]$id) #bacterial IDs
correct <- c(correct, results[[i]]$correct.id)
}

pvalue <- p_estimate(log.probs, correct) #estimate p-values using decoy results
fdr <- p.adjust(pvalue, method="fdr") #remove decoy results

sum(fdr[correct>=0]<=0.05) #180
sum(correct[fdr<=0.05]==0)/sum(fdr[correct>=0]<=0.05) #2.22%

