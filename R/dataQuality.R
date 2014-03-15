# functions comparing the data sets and measuring quality of generated data


dataSimilarity <- function(data1, data2, dropDiscrete=NA) {
  
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  
  
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  data1Num <- data1[,!discrete,drop=F]
  data2Num <- data2[,!discrete,drop=F]
  
  s1<- matrix(0,nrow=4,ncol=ncol(data1Num))
  s2<- matrix(0,nrow=4,ncol=ncol(data1Num))
  s1Norm<- matrix(0,nrow=4,ncol=ncol(data1Num))
  s2Norm<- matrix(0,nrow=4,ncol=ncol(data1Num))
  ks <- c() # store p-values of Kolmogorov-Smirnov test on drawing from the same distribution, only for numeric attributes
  
  
  if (ncol(data1Num) > 0) {
    
    s1[1,] <- colStats(data1Num, mean, na.rm=TRUE)
    s1[2,] <- colStats(data1Num, sd, na.rm=TRUE)
    s1[3,] <- colStats(data1Num, skewness, na.rm=TRUE)
    s1[4,] <- colStats(data1Num, kurtosis, na.rm=TRUE)
    
    s2[1,] <- colStats(data2Num, mean, na.rm=TRUE)
    s2[2,] <- colStats(data2Num, sd, na.rm=TRUE)
    s2[3,] <- colStats(data2Num, skewness, na.rm=TRUE)
    s2[4,] <- colStats(data2Num, kurtosis, na.rm=TRUE)
    
    # normalize data to [0,1] 
    data1NumNorm <- normalize01(data1Num)
    data2NumNorm <- normalize01(data2Num)
    
    s1Norm[1,] <- colStats(data1NumNorm, mean, na.rm=TRUE)
    s1Norm[2,] <- colStats(data1NumNorm, sd, na.rm=TRUE)
    s1Norm[3,] <- colStats(data1NumNorm, skewness, na.rm=TRUE)
    s1Norm[4,] <- colStats(data1NumNorm, kurtosis, na.rm=TRUE)
    
    s2Norm[1,] <- colStats(data2NumNorm, mean, na.rm=TRUE)
    s2Norm[2,] <- colStats(data2NumNorm, sd, na.rm=TRUE)
    s2Norm[3,] <- colStats(data2NumNorm, skewness, na.rm=TRUE)
    s2Norm[4,] <- colStats(data2NumNorm, kurtosis, na.rm=TRUE)
    
    for (i in 1:ncol(data1Num)) { # KS test 
      tst <- ks.test(data1Num[,i], data2Num[,i])
      ks[i] <- tst$p.value
    }
    names(ks) <- names(data1Num)
    
  }
  # difference between normalized statistics
  ds <- s1Norm - s2Norm
  
  colnames(s1) <- colnames(s2) <- colnames(s1Norm) <- colnames(s2Norm) <- colnames(ds) <- names(data1Num)
  rownames(s1) <- rownames(s2) <- rownames(s1Norm)<- rownames(s2Norm) <- rownames(ds) <- c("mean","stdev","skewness","kurtosis")
  
  # nominal attributes: compare frequencies and Hellinger distance
  f1 <- list()
  f2 <- list()
  df <- list()
  hel <- c() # store Hellinger distance between two discrete distributions
  discreteIdx<-which(discrete)
  if (! is.na(dropDiscrete)) {
    dropIdx <- match(dropDiscrete, names(data1))
    discreteIdx <- discreteIdx[-match(dropIdx, discreteIdx, nomatch=length(discreteIdx)+1)]
  }
  for (i in seq(along=discreteIdx)) {
    f1[[i]] <- table(data1[,discreteIdx[i]])/nrow(data1)
    names(f1)[i] <-names(data1)[discreteIdx[i]]
    f2[[i]] <- table(data2[,discreteIdx[i]])/nrow(data2)
    names(f2)[i] <-names(data2)[discreteIdx[i]]
    df[[i]] <- f1[[i]] - f2[[i]]
    names(df)[i] <-names(data1)[discreteIdx[i]]
    hel[i] <- hellinger(f1[[i]], f2[[i]])
  }
  if (length(hel) > 0) 
      names(hel)<- names(data1)[discreteIdx]
  
  list(equalInstances=equalInstDaisy(data1,data2), stats1num=s1, stats2num=s2, stats1numNorm=s1Norm, stats2numNorm=s2Norm,
       KSpvalue=ks, freq1=f1, freq2=f2, dfreq=df, dstatsNorm=ds, hellingerDist=hel)
}

# equalInst<-function(data1, data2) {
#   eq <- 0
#   for (i in 1:nrow(data2)) {
#     row2 <- data2[i,]
#     eq <- eq + sum(apply(data1, 1, almostEqual, inst2=row2, tolerance=1e-5))
#   }
#   eq
# }
# 
equalInst<-function(data1, data2) {
  eq <- 0
  for (i in 1:nrow(data1)) 
    for (j in 1:nrow(data2))    
      eq <- eq + as.integer(almostEqual(data1[i,], data2[j,],tolerance=1e-5))
  eq
}

equalInstDaisy <- function(data1, data2, tolerance=1e-5) {
  dta <- rbind(data1, data2)
  sim <- as.matrix(daisy(dta, metric="gower"))
  eq <- 0
  for (i in (nrow(data1)+1):nrow(sim)) {
    if (length(which(sim[i,]< tolerance)) >1)
      eq <- eq+1
  }
  eq
}

meanMinDistance<-function(data1, data2) {
  dist <- vector(mode="numeric",length=nrow(data2))
  for (i in 1:nrow(data2)) {
    row2 <- as.numeric(data2[i,])
    dist[i] <- min(apply(data1, 1, diff, row2))
  }
  
  list(mean=mean(dist),min=min(dist),dist=dist)
}

diff<-function(inst1, inst2){
  sqrt(sum((inst1-inst2)^2,na.rm=T))
}

almostEqual <- function(inst1, inst2, tolerance=1e-5) {
  all(abs(as.numeric(inst1)-as.numeric(inst2)) < tolerance, na.rm=TRUE)
}

# compare adjusted Rand index of two cluserings obtained by:
# 1.  using data1 to get clusters, assign data2 to these clusters
# 2. using data2 to get clusters, assign data1 to these clusters
# The method uses pam clustering (partitioning around medoids) with 
# optimal k computed for each data1 and data2 separately with average silhouette width method (using pamk)
# 
# ARI is computed on merged data1 and data2 instances 
ariCompare <- function(data1, data2) {
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  #distance matrix
  dist1 <- daisy(data1, metric="gower")
  dist2 <- daisy(data2, metric="gower")
  
  pamkBest1 <- pamk(dist1) # number of clusters estimated by optimum average silhouette width
  pamkBest2 <- pamk(dist2) # number of clusters estimated by optimum average silhouette width
  
  clu1 <- pam(dist1, k=max(2,pamkBest1$nc)) #get the clustering
  clu2 <- pam(dist2, k=max(2,pamkBest2$nc)) #get the clustering
  
  medoids1 <- data1[clu1$medoids,] 
  medoids2 <- data2[clu2$medoids,] 
  
  # compute distances of data2 to medoids1 and data1 to medoids2
  dta1 <- rbind(medoids1,data2)
  dta2 <- rbind(medoids2,data1)
  dst1 <- as.matrix(daisy(dta1,metric="gower"))
  dst2 <- as.matrix(daisy(dta2,metric="gower"))
  # find nearest medoid for each data2 instance, this is its cluster index in 1st clustering
  closestMed1 <-apply(dst1[-c(1:nrow(medoids1)),1:nrow(medoids1)],1,which.min)
  # find nearest medoid for each data1 instance, this is its cluster index in 2nd clustering
  closestMed2 <-apply(dst2[-c(1:nrow(medoids2)),1:nrow(medoids2)],1,which.min)
  
  # put the cluster assignments together
  clustering1 <- c(clu1$clustering, closestMed1)
  clustering2 <- c(closestMed2, clu2$clustering)
  
  adjustedRandIndex(clustering1, clustering2)
}

performanceCompare <- function(data1, data2, formula, model="rf", stat="accuracy", ...) {
  
  informula <- as.formula(formula)
  trms <- terms(as.formula(informula),data=data1)
  response <- all.names(attr(trms,"variables"))[1+attr(trms,"response")]
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  # shuffle the data to allow several different runs
  d1 <- data1[sample(nrow(data1),nrow(data1)),]
  d2 <- data2[sample(nrow(data2),nrow(data2)),]
  
  # split into halves based on startified cross-validation
  cv1 <- cvGenStratified(d1[[response]], k=2)
  d1a <- d1[cv1==1,]
  d1b <- d1[cv1==2,]
  
  cv2 <- cvGenStratified(d2[[response]], k=2)
  d2a <- d2[cv2==1,]
  d2b <- d2[cv2==2,]
  
  # build models
  m1a <- CoreModel(informula, d1a, model=model, ...)
  m1b <- CoreModel(informula, d1b, model=model, ...)
  
  m2a <- CoreModel(informula, d2a, model=model, ...)
  m2b <- CoreModel(informula, d2b, model=model, ...)
  
  # test performance
  # model m1a
  pred.m1ad1b <- predict(m1a, d1b)
  acc.m1ad1b <- modelEval(m1a, d1b[[response]], pred.m1ad1b$class, pred.m1ad1b$probabilities)  
  
  pred.m1ad2 <- predict(m1a, d2)
  acc.m1ad2 <- modelEval(m1a, d2[[response]], pred.m1ad2$class, pred.m1ad2$probabilities)  
  
  #model m1b
  pred.m1bd1a <- predict(m1b, d1a)
  acc.m1bd1a <- modelEval(m1b, d1a[[response]], pred.m1bd1a$class, pred.m1bd1a$probabilities)  
  
  pred.m1bd2 <- predict(m1b, d2)
  acc.m1bd2 <- modelEval(m1b, d2[[response]], pred.m1bd2$class, pred.m1bd2$probabilities)  
  
  # model m2a
  pred.m2ad2b <- predict(m2a, d2b)
  acc.m2ad2b <- modelEval(m2a, d2b[[response]], pred.m2ad2b$class, pred.m2ad2b$probabilities)  
  
  pred.m2ad1 <- predict(m2a, d1)
  acc.m2ad1 <- modelEval(m2a, d1[[response]], pred.m2ad1$class, pred.m2ad1$probabilities)  
  
  #model m2b
  pred.m2bd2a <- predict(m2b, d2a)
  acc.m2bd2a <- modelEval(m2b, d2a[[response]], pred.m2bd2a$class, pred.m2bd2a$probabilities)  
  
  pred.m2bd1 <- predict(m2b, d1)
  acc.m2bd1 <- modelEval(m2b, d1[[response]], pred.m2bd1$class, pred.m2bd1$probabilities)  
  
  
  # averages for 2-fold cross-validation
  acc.m1d1 <- (acc.m1ad1b[[stat]] + acc.m1bd1a[[stat]])/2.0
  acc.m1d2 <- (acc.m1ad2[[stat]] + acc.m1bd2[[stat]])/2.0
  acc.m2d1 <- (acc.m2ad1[[stat]] + acc.m2bd1[[stat]])/2.0
  acc.m2d2 <- (acc.m2ad2b[[stat]] + acc.m2bd2a[[stat]])/2.0
  diff.m1 <- acc.m1d1 -acc.m1d2
  diff.m2 <- acc.m2d2 - acc.m2d1
  
  list(diff.m1=diff.m1, diff.m2=diff.m2, perf.m1d1=acc.m1d1, perf.m1d2=acc.m1d2, perf.m2d1=acc.m2d1, perf.m2d2=acc.m2d2 )
  
}


performanceCompareSingle <- function(data1, data2, formula, model="rf", stat="accuracy", ...) {
  
  informula <- as.formula(formula)
  trms <- terms(as.formula(informula),data=data1)
  response <- all.names(attr(trms,"variables"))[1+attr(trms,"response")]
  
  # check if data is of compatibile types
  if(ncol(data1)!=ncol(data2))
    stop("Only data with equal number of columns can be compared.")
  discrete <- c()
  noVal <- c()
  for (i in 1:ncol(data1)) {
    if (is.factor(data1[[i]]) != is.factor(data2[[i]]))
      stop("Only data with equal types of columns can be compared.", names(data1)[i], " ", names(data2)[i])
    discrete[i] <- is.factor(data1[[i]])
    if (discrete[i]) {
      if (! all(levels(data1[[i]]) == levels(data2[[i]])))
        stop("Only data with equal types of columns can be compared.", levels(data1[[i]]), " ", levels(data2[[i]]))
      noVal[i] <- length(levels(data1[[i]]))
    }
    else
      noVal[i] <- NA
  }
  
  # shuffle the data
  d1 <- data1[sample(nrow(data1),nrow(data1)),]
  d2 <- data2[sample(nrow(data2),nrow(data2)),]
  
   
  # build models
  m1 <- CoreModel(informula, d1, model=model, ...)
  m2 <- CoreModel(informula, d2, model=model, ...)

  # test performance
  # model m1
  pred.m1d2 <- predict(m1, d2)
  acc.m1d2 <- modelEval(m1, d2[[response]], pred.m1d2$class, pred.m1d2$probabilities)  
  
  # model m2a
  pred.m2d1 <- predict(m2, d1)
  acc.m2d1 <- modelEval(m2, d1[[response]], pred.m2d1$class, pred.m2d1$probabilities)  
  
  
  # difference in performance
  diff <- acc.m1d2[[stat]] - acc.m2d1[[stat]]
  
  list(m1=acc.m1d2[[stat]], m2=acc.m2d1[[stat]], diff=diff)
}



 