#Author: So Young Ryu
#Date: 03/25/2021
#Functions for BacID Score and BacID Prob
#Ryu, S., Likelihood-based Bacterial Identification approach for bimicrobial mass spectrometry data (Revised and Submitted to Annal of Applied Statistics)

#BacID Prob Function
#poly: an observed bimicrobial mass spectrum
#ref: a reference data which contains spectra from both real species and decoy
#true.id: Please enter true.id for spectra from known bacterial mixtures
#mass.tolerance (in ppm): 2000 by default
#mz.sd0: mz.sd0 where peak m/z standard deviation is mz.sd0 times a reference m/z
#int.sd: peak intensity standard deviation
#ini.gamma: an initial gamma value, 0.5 by default
#ini.alpha: a vector of alphas (length 2). If not given, an initial alphas will be estimated.
#esp: a EM algorithm convergence criteria
#max.iter: the maximum number of iterations for EM algorithm
#Return the followings:
#id.paired: Top paired IDs for the model (will be used for BacID Score)
#id: two bacterial identifications (IDs for BacID Prob)
#correct.id: If TRUE.ID is given, this is an indicator variable that IDs is correct IDs(1), incorrect IDs(0) from original database or decoy IDs (01). If TRUE.ID is not given, 0 is for original IDs and -1 is for decoy IDs. 
#total.iter: the number of EM iterations
#log.probs: log of posterior probabilties for identified bacteria

matching_bi <- function(poly, ref,   true.id=NULL, mass.tolerance=2000,
						mz.sd0=0.00055, int.sd=0.6, ini.gamma=0.50,
						ini.alpha=NULL, esp=0.001, max.iter=50){
	m=2
	db <- combn(names(ref),m)
	loglik.list = prob.list =  array(c(NA), dim=ncol(db))
		for (i2 in 1:ncol(db)){ 
		index <- which(names(ref) %in% db[,i2])
		refs <- ref[index]
		monos.ini = monos = list()
		aligned.all <- poly
			for (i3 in 1:length(refs)){
			monos.ini[[i3]] <- data.frame(mz=refs[[i3]][,1], int=refs[[i3]][,2])		
			common.peaks <- find.common.peaks(poly, monos.ini[[i3]], tolerance=mass.tolerance)
			aligned <- align.peaks(common.peaks); aligned <- aligned[!is.na(aligned[,1]),]
			mono <- data.frame(mz=aligned[,3], int=aligned[,4])
			monos[[i3]] <- mono	
			aligned.all <- cbind(aligned.all, mono)
			}
		res <- fit_model_C2(poly, monos, mass.tolerance=mass.tolerance,
				mz.sd0=mz.sd0, int.sd=int.sd, ini.gamma=ini.gamma,
				ini.alpha=ini.alpha, esp=esp, max.iter=max.iter)
		loglik.list[i2] <- res$loglik
		prob.list[i2] <- exp(loglik.list[i2])
		if (prob.list[i2]>1 || prob.list[i2]<0){ cat(i2, prob.list[i2], "\n")}
		}
	  if (sum(prob.list)>0){
		prob.list=prob.list/sum(prob.list)}
		
		ids.paired <- db[,which.max(loglik.list)]
		
		probs.ref <- rep(0, length(ref))
		for (i3 in 1:ncol(db)){
			for (i4 in 1:length(ref)){
				if (names(ref)[i4] %in% db[,i3]){
				probs.ref[i4]=probs.ref[i4]+prob.list[i3]	
				}
			}
		}
		order <- sort.int(probs.ref, decreasing=TRUE, index.return=TRUE)$ix
	#numerical errors occur sometimes
	if (probs.ref[order[1]] > (1-1e-10)){probs.ref[order[1]]=1}
	if (probs.ref[order[2]] > (1-1e-10)){probs.ref[order[2]]=1}
	id1 <- names(ref)[order[1]]; log.prob1 <- log(probs.ref[order[1]])
	id2 <- names(ref)[order[2]]; log.prob2 <- log(probs.ref[order[2]]) 
	correct.id1 <- ifelse(id1 %in% true.id, 1, ifelse(substr(id1, 1,5)=='decoy',-1,0))	
	correct.id2 <- ifelse(id2 %in% true.id, 1, ifelse(substr(id2, 1,5)=='decoy',-1,0))	

	return(list( ids.paired = ids.paired, id=c(id1, id2), correct.id=c(correct.id1, correct.id2), total.iter=res$total.iter, log.probs=c(log.prob1, log.prob2)))
}



#BacID Score Function
#poly: an observed bimicrobial mass spectrum
#ref: a reference data which contains spectra from both real species and decoy
#true.id: Please enter true.id for spectra from known bacterial mixtures
#m: The total number of resampling
#p: a resampling proportion
#mass.tolerance (in ppm): 2000 by default
#mz.sd0: mz.sd0 where peak m/z standard deviation is mz.sd0 times a reference m/z
#int.sd: peak intensity standard deviation
#ini.gamma: an initial gamma value, 0.5 by default
#ini.alpha: a vector of alphas (length 2). If not given, an initial alphas will be estimated.
#esp: a EM algorithm convergence criteria
#max.iter: the maximum number of iterations for EM algorithm	
#Return the followings:
#id: two bacterial identifications
#correct.id: If TRUE.ID is given, an indicator variable saying that IDs is correct IDs(1), incorrect IDs(0) from original database or decoy IDs (01). If TRUE.ID is not given, 0 for original IDs and -1 for decoy IDs. 
#scores: BacID scores for two identified bacteria

conf_score <- function(poly, ref, true.id, m=100, p=0.70, mass.tolerance=2000,
				mz.sd0=0.00055, int.sd=0.6, ini.gamma=0.50,
				ini.alpha=NULL, esp=0.001, max.iter=50){
	n=ceiling(nrow(poly)*p)
	score <- array(c(0), dim=length(ref))
	names(score) <- names(ref)
	for (i in 1:m){
		index=sort(sample(1:nrow(poly))[1:n])
		poly.j <- poly[index,]	
		id.j <- matching_bi(poly.j, ref, true.id, mass.tolerance=mass.tolerance,
						mz.sd0=mz.sd0, int.sd=int.sd, ini.gamma=ini.gamma,
						ini.alpha=ini.alpha, esp=esp, max.iter=max.iter)$ids.paired
		bingo.index <- which(names(ref) %in% id.j)
		score[bingo.index]=score[bingo.index]+1
	}
	score=score/m
	
	order <- sort.int(score, decreasing=TRUE, index.return=TRUE)$ix
	id1 <- names(ref)[order[1]]; score1 <- score[order[1]]
	id2 <- names(ref)[order[2]]; score2 <- score[order[2]] 
	correct.id1 <- ifelse(id1 %in% true.id, 1, ifelse(substr(id1, 1,5)=='decoy',-1,0))	
	correct.id2 <- ifelse(id2 %in% true.id, 1, ifelse(substr(id2, 1,5)=='decoy',-1,0))		
	return(list(ids=c(id1, id2), correct.id=c(correct.id1,correct.id2), scores=c(score1, score2)))
}



#Fit Equation (2) using EM algorithm.
#poly: an observed bimicrobial mass spectrum
#monos: reference spectra of a size two
#mass.tolerance (in ppm): 2000 by default
#mz.sd0: mz.sd0 where peak m/z standard deviation is mz.sd0 times a reference m/z
#int.sd: peak intensity standard deviation
#ini.gamma: an initial gamma value, 0.5 by default
#ini.alpha: a vector of alphas (length 2). If not given, an initial alphas will be estimated.
#Epsilon: a EM algorithm convergence criteria
#max.iter: the maximum number of iterations for EM algorithm
#Return the followings:
#param: Parameters (alpha1, alpha2, gamma)
#loglik: observed data log likelihood
#loglik.trace: a trace of observed data log likelihood
#e.q: probabilities that a peak belongs to the first bacteria
#total.iter: Total number of iterations  

fit_model_C2 <- function(poly, monos, mass.tolerance=2000,
				mz.sd0=0.00055, int.sd=0.6, ini.gamma=0.50,
				ini.alpha=NULL, esp=0.001, max.iter=50){ 
					
conv.fun <- function(old, new, esp) abs(new-old) < esp*(1+abs(new))
conv <- function(l.old, l.new, param.old, param.new, status, esp){
	conv.l <- conv.fun(l.old, l.new, esp)
	conv.alpha1 <- conv.fun(param.old$alpha[1], param.new$alpha[1], esp)
	conv.alpha2 <- conv.fun(param.old$alpha[2], param.new$alpha[2], esp)
	if (sum(status==3)>0){
	conv.gamma <- conv.fun(param.old$gamma, param.new$gamma, esp)
	}else{conv.gamma=TRUE} #no need to compute gamma
	if (conv.l && conv.alpha1 && conv.alpha2 && conv.gamma){return(TRUE)
	#if (conv.l){return(TRUE)
	}else{return(FALSE)}
}
e_step <- function(d.mz, d.int, status, param,  int.sd){
	e.q = rep(NA, length(status))
	e.q[status==1]=1
	e.q[status==2]=0
	if (sum(status==3)>0){
	num = param$gamma*f(d.mz[[1]][status==3], d.int[[1]][status==3], param$alpha[1], int.sd)
	den = param$gamma*f(d.mz[[1]][status==3], d.int[[1]][status==3], param$alpha[1], int.sd)+(1-param$gamma)*f(d.mz[[2]][status==3], d.int[[2]][status==3], param$alpha[2], int.sd)
	e.q[status==3]=(num/den)
	}
	return(e.q)
}

f <- function(d.mz, d.int, alpha, int.sd){
	res=dnorm(d.mz, mean=0, sd=1)*dnorm(d.int, mean=alpha/int.sd, sd=1)
	return(res)
}
m_step <- function(d.mz, d.int, status, e.q, int.sd){
	gamma <- sum(e.q[status==3])/sum(status==3)	
	alpha <- c()
	alpha[1] <- ifelse((sum(e.q[status==3])+sum(status==1))>0,int.sd*(sum(e.q[status==3]*d.int[[1]][status==3])+sum(d.int[[1]][status==1]))/(sum(e.q[status==3])+sum(status==1)),0)
	alpha[2] <- ifelse((sum(1-e.q[status==3])+sum(status==2))>0, int.sd*(sum((1-e.q[status==3])*d.int[[2]][status==3])+sum(d.int[[2]][status==2]))/(sum(1-e.q[status==3])+sum(status==2)),0)
	return(list(gamma=gamma, alpha=alpha))
}
loglik <- function(d.mz, d.int, status, param, r.int, r.mz, int.sd){
	log.lik = -(log(r.mz)+log(r.int))*sum(status==0)+
	sum(log(param$gamma*f(d.mz[[1]][status==3], d.int[[1]][status==3], param$alpha[1], int.sd)+(1-param$gamma)*f(d.mz[[2]][status==3], d.int[[2]][status==3], param$alpha[2], int.sd)))+
	sum(log(f(d.mz[[1]][status==1], d.int[[1]][status==1], param$alpha[1], int.sd)))+
	sum(log(f(d.mz[[2]][status==2], d.int[[2]][status==2], param$alpha[2], int.sd)))
	return(log.lik)
}
				
###############
#Main Function#
###############
min.count=3
r.mz <- (12000-4000)/(4000*mz.sd0)
r.int <- max(poly$int)/int.sd
if (sum(monos[[1]]$int>0)+sum(monos[[2]]$int>0) < min.count){
	l = rep(1/r.mz/r.int, nrow(poly))
	log.l <- sum(log(l))
	return(list(loglik=log.l))
}	
status <- ifelse(monos[[1]]$int==0 & monos[[2]]$int==0, 0, #no peak is matched
				ifelse(monos[[1]]$int!=0 & monos[[2]]$int==0, 1, #ref1 matched
				ifelse(monos[[1]]$int==0 & monos[[2]]$int!=0, 2, #ref2 matched
				3))) #peaks from both refs are matched
if (length(ini.alpha)==0){
	ini.alpha[1] = ifelse(sum(status==1)>0, median((poly$int-monos[[1]]$int)[status==1]), 0)
	ini.alpha[2] = ifelse(sum(status==2)>0, median((poly$int-monos[[2]]$int)[status==2]),0)	
}
d.mz <- list(); d.int <- list()
for (k in 1:2){
	d.mz[[k]]=(poly$mz-monos[[k]]$mz)/(monos[[k]]$mz*mz.sd0)
	d.int[[k]]=(poly$int-monos[[k]]$int)/int.sd
}
param <- list(); 
param$alpha <- ini.alpha
param$gamma <- ini.gamma				
loglik.trace <- NULL
to.run <- TRUE
j <- 0
	while(to.run){
		j <- j+1				
		e.q <- e_step(d.mz, d.int, status, param,  int.sd)
		if (j>1) {param.old <- param}
		param <- m_step(d.mz, d.int, status, e.q,  int.sd)
		if (j>1) {l.old <- l.new}
		l.new <- loglik(d.mz, d.int, status, param, r.int, r.mz, int.sd)
		loglik.trace <- c(loglik.trace, l.new)
		if (j==1) {to.run=TRUE
		}else{
		if (j>=max.iter) {to.run <- FALSE
		}else {to.run <- !conv(l.old, l.new, param.old, param, status, esp)}
	} 
}
return(list(param=param, loglik=l.new,  loglik.trace=loglik.trace,e.q=e.q, total.iter=j))	
}				


									
#normalize intensity
#data: mass spectra
#log.transform: log-transformation of peak intensities, TRUE by default
#intensity: 90th percentile intensity
#constant: add a small constant to avoid log(negative #)
#Return a transformed data (mass spectra with m/z and intensity columns)
normalize.intensity = function(data, log.transform=TRUE, intensity = 100, constant = 1) {
  new.data = data
  new.data[, 2] = new.data[, 2] / quantile(new.data[, 2], prob=0.90) * intensity
  if (log.transform==TRUE){
    new.data[,2]=log(new.data[,2]+constant)}
  new.data
}



#estimate p vlaues using decoy results
#probs: posterior probabilities or scores
#org: a vector indicating that it is from an original database or decoy database. 1 for original -1 for decoy
p_estimate <- function(probs, org){	
p <- c()
for (i1 in 1:length(probs)){
p[i1] = sum(org[probs>=probs[i1]]==-1)/sum(org>=0)
}
return(p)	
}



#Codes from Yang et al. (2018) Direct MALDI-TOF MS Identification of Bacterial Mixtures, Analytical Chemistry. 
#Find common peaks between two spectra 
#peaks1: paired m/z and intensity values from the first spectrum
#peaks2: paired m/z and intensity values from the second spectrum
#tolerance: mass tolerance 
find.common.peaks = function(peaks1, peaks2, tolerance) {
  if(class(peaks1) == 'data.frame')
    peaks1 = as.matrix(peaks1)
  if(class(peaks2) == 'data.frame')
    peaks2 = as.matrix(peaks2)
  
  indexs.1 = order(peaks1[, 2], decreasing = TRUE)
  
  is.used = rep(FALSE, nrow(peaks2))
  indexs.2 = sapply(indexs.1, function(i) {
    mz1 = peaks1[i, 1]
    candidates = which(!is.used & peaks2[, 1] >= mz1 * (1 - tolerance * 1e-6) & peaks2[, 1] <= mz1 * (1 + tolerance * 1e-6))
    j = candidates[which.min(abs(peaks2[candidates, 1] - mz1))]
    if(length(j) == 0)
      return(NA)
    else {
      is.used[j] <<- TRUE
      return(j)
    }
  })
  
  common.indexs = which(!is.na(indexs.2))
  if(length(common.indexs) == 0)
    list(
      common.peaks1 = peaks1[c(), ],
      common.peaks2 = peaks2[c(), ],
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
  else if(length(common.indexs) == 1)
    list(
      common.peaks1 = matrix(peaks1[indexs.1[common.indexs], ], nrow = 1),
      common.peaks2 = matrix(peaks2[indexs.2[common.indexs], ], nrow = 1),
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
  else
    list(
      common.peaks1 = peaks1[indexs.1[common.indexs], ],
      common.peaks2 = peaks2[indexs.2[common.indexs], ],
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
}



#Codes from Yang et al. (2018) Direct MALDI-TOF MS Identification of Bacterial Mixtures, Analytical Chemistry.  
#Aligned Peaks an object, "common.peaks" resulted from a function "find.common.peaks"
#common.peaks: an output of a function of "find.common.peaks"
align.peaks = function(common.peaks) {
  peaks1.common = cbind(common.peaks[[1]][, 1], common.peaks[[1]][, 2])
  peaks2.common = cbind(common.peaks[[2]][, 1], common.peaks[[2]][, 2])
  peaks1 = cbind(common.peaks[[3]][, 1], common.peaks[[3]][, 2])
  peaks2 = cbind(common.peaks[[4]][, 1], common.peaks[[4]][, 2])
  if(nrow(peaks1.common) > 0) {
    get.noncommom.peaks.index = function(peaks, peaks.common) {
      sapply(1:nrow(peaks), function(i) {
        idx = which(peaks.common[, 1] == peaks[i, 1] & peaks.common[, 2] == peaks[i, 2])
        if(length(idx > 0)) {
          peaks.common[idx[1], 1] <<- NA
          FALSE
        }
        else {
          TRUE
        }
      })
    }
    
    mz2.noncommon = peaks2[get.noncommom.peaks.index(peaks2, peaks2.common), 1]
    mz1.noncommon = peaks1[get.noncommom.peaks.index(peaks1, peaks1.common), 1]
    peaks1.new = if(length(mz2.noncommon) > 0) rbind(peaks1, cbind(mz2.noncommon, 0)) else peaks1
    peaks2.new = if(length(mz1.noncommon) > 0) rbind(peaks2, cbind(mz1.noncommon, 0)) else peaks2
  }
  else {
    peaks1.new = rbind(peaks1, cbind(peaks2[, 1], 0))
    peaks2.new = rbind(peaks2, cbind(peaks1[, 1], 0))
  }
  peaks1.new = peaks1.new[order(peaks1.new[, 1]), ]
  peaks2.new = peaks2.new[order(peaks2.new[, 1]), ]
  peaks1.new[which(peaks1.new[, 2] == 0), 1] = NA
  peaks2.new[which(peaks2.new[, 2] == 0), 1] = NA
  colnames(peaks1.new) = c('mz1', 'int1')
  colnames(peaks2.new) = c('mz2', 'int2')
  return(cbind(peaks1.new, peaks2.new))
}



#mass recalibration when there are systemic errors in m/z measurement per species
#ref: a reference database
#mz.shift: a vector of m/z value correction value (a positive value when a reference spectrum has a higher m/z value compared to an observed value)
#Return a mass recalibrated reference database
mass.recal <- function(ref, mz.shift){
  for (i in 1:length(mz.shift)){
  ref[[i]][,1] = ref[[i]][,1]+mz.shift[i]
  }
  return(ref)
}





		
		
				
				