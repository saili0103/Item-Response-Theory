
DIF_logit <- function(theta, item, diag){
  #item = halluc.mat[,i];theta = theta.est; diag = group0
  #if(spline.fit == T){
   # theta <- spline(1:length(theta),theta)$y
    #plot(1:length(theta),theta)
    #lines(spline(1:length(theta),theta))
  #}
  not.na <- as.vector(which(!is.na(item)))
  item <- item[not.na]
  theta <- theta[not.na]
  item <- as.factor(item)
  diag <- as.factor(diag[not.na])
  diag <- droplevels(diag)
  if(nlevels(diag)< 2){
    return(NA)
  }
  data1 <- data.frame(item = item,theta = theta, diag = diag)

  if(nlevels(item) > 2){
    fit1 <- tryCatch(polr(item ~ theta + diag + theta*diag, data = data1), error = function(e) NULL)
    fit2 <- tryCatch(polr(item ~ theta, data = data1), error = function(e) NULL)
    if(is.null(fit1) | is.null(fit2)){
      return(NA)
    }
    anov.re = anova(fit1, fit2)
    deltaG <- anov.re$LR[2]
    dfG <- anov.re$`   Df`[2]
  }else{
    logit.fit <-tryCatch(glm(item ~ theta + diag + diag * theta, data = data1, family= 'binomial'), error = function(e) NULL)
    if(is.null(logit.fit)){
      return(NA)
    }
    anov.re <- anova(logit.fit)
    deltaG <- anov.re$Dev[3] + anov.re$Dev[4]
    dfG <-anov.re$Df[3] + anov.re$Df[4]
    #statistic :P332
  }
  N = length(item)
  power.chisq <- pwr.chisq.test(deltaG/N, N= N, df = dfG, sig.level = 0.05)$power
  #critical.value <- c(qchisq(0.95, dfG[1],lower.tail = T), qchisq(0.95, dfG[2], lower.tail = T))
  list(stat = deltaG, df = dfG, power = power.chisq, N = N)
}

DIF_MH <- function(theta, item, diag){
  #item = halluc.mat[,i];theta = theta.est; diag = group0

  not.na <- as.vector(which(!is.na(item)))
  item <- item[not.na]
  theta <- theta[not.na]
  item <- as.factor(item)
  diag <- as.factor(diag[not.na])
  diag <- droplevels(diag)
  if(nlevels(diag)< 2){
    return(NA)
  }
  N <- length(item)
  data1 <- data.frame(item = item,theta = theta, diag = diag)
  nlevels <- nlevels(item)
  qth <- quantile(theta, seq(0, 1, by = .1), na.rm = TRUE)
  while(any((qth[2:11] - qth[1:10]) < 10^(-10))){
    theta <- theta + rnorm(length(theta), mean = 0, sd = 10^(-4))
    qth <- quantile(theta, seq(0, 1, by = .1), na.rm = TRUE)
  }
  qth[1] <- qth[1] - 1e-06
  qth[11] <- qth[11] + 1e-06

  groups.ind <- cut(theta, breaks = qth)
  no.group <- nlevels(groups.ind)
  levels(groups.ind) <- as.character(1: no.group)
 # Nj <- table(groups.ind)
  array.vec <- NULL
  for(k in 1: no.group){
    data.k <- which(groups.ind == k)
    array.vec <- c(array.vec, as.vector(table(item[data.k], diag[data.k])))
  }
  array.MH <-
    array(array.vec,
          dim = c(nlevels(item), 3, 10),
          dimnames = list( Response = levels(item), Diag = c("DSM4295", "DSM4V71", "DSM4296"),
            theta.group = as.character(1:10)))
  ## Classical Mantel-Haenszel test
  mh.re <- tryCatch(mantelhaen.test(array.MH),error = function(e) NULL)
  if(is.null(mh.re)){
    return(NA)
  }
  power.chisq <- pwr.chisq.test(mh.re$stat/N, N= N, df = mh.re$par, sig.level = 0.05)$power
  list(stat = mh.re$stat, N = N, df = mh.re$par, power = power.chisq)
}