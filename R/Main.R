#' Estimate the value of theta_0 using MLE with ridge penalty
#'
#' @param nodesLabel The vector of latent state of all genes.
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param lambda The tuning parameter for ridge.
#' @param threshold The threshold for convergence.
#' @param maxite The maximum number of iteration.
#' @return The estimated value for theta_0:\code{beta0_new}
#'

Get_theta_0<-function(nodesLabel,nodesNeighbor,w,W,lambda=0.1,threshold=1e-4,maxite=100){
  if(!is.null(W)){
    print("Weighted matrix!")
    x <- cbind(x0=1,x1=sapply(nodesNeighbor,function(neiIdx)  sum(W[nodesLabel[neiIdx]== 1,])),
               x2=sapply(nodesNeighbor,function(neiIdx) -sum(W[nodesLabel[neiIdx]==-1,])))
  }else{
    x <- cbind(x0=1,x1=mapply(function(neiIdx,wi)  sum(wi+w[neiIdx][nodesLabel[neiIdx]== 1]),nodesNeighbor,w),
               x2=mapply(function(neiIdx,wi) -sum(wi+w[neiIdx][nodesLabel[neiIdx]==-1]),nodesNeighbor,w))
  }
  beta0 <-c(0,0,0)
  y <- (nodesLabel+1)/2
  np <- ncol(x)
  converg <- rep(F,np)
  ite <- 0
  while(!all(converg) && ite<maxite) {
    ite <- ite+1
    linearpred <- as.numeric(x %*% beta0)
    prob <- 1/(1+exp(- linearpred))
    beta0_new <- beta0 + solve(t(x) %*% diag(prob*(1-prob)) %*% x  + diag(2*lambda,np)) %*% ( t(x)%*% (y-prob) - 2*lambda*beta0 )
    converg <- abs(beta0_new-beta0)/(0.1+abs(beta0_new)) < threshold
    beta0 <- beta0_new
    # browser()
    cat("beta0:", beta0_new,"\n")
  }
  if(!all(converg)){print(paste0("Logistic regression for theta_0 did not converge after ",maxite," iterations!"))}else{
    print(paste0("Logistic regression for theta_0 converges after ",ite," iterations!"))
  }

  beta0_new<-as.numeric(beta0_new)
  names(beta0_new)<-c("h1","tau_11","tau_00")
  return(beta0_new)
}

#' Estimate for estimating data parameters theta_1
#'
#' @param nodesLabel Latent state of all genes.
#' @param Y The vector of DNM counts for all genes.
#' @param N The parameter in the prior for network.
#' @param mu The vector of mutability for all genes.
#' @param beta_1_init The initiate value for Newton's method.
#' @return The estimated value for theta_1:\code{beta_1_new}
#'
Get_theta_1_Gamma<-function(nodesLabel,Y,N,mu,beta_1_init){
  y <- (nodesLabel+1)/2

  Y_sub<-Y[which(y==1)]
  y_sub<-y[which(y==1)]
  mu_sub<-mu[which(y==1)]

  beta_1 <- 99999
  beta_1_new <- beta_1_init

  #update beta
  beta_1 <- beta_1_new
  beta_1_new <- log(sum(Y_sub) / sum(2*N*mu_sub))
  print(paste0("Estimation of theta_1 (gamma) is ",exp(beta_1_new)))
  return(beta_1_new)
}



#'Calculate the conditional probability: P(Si|S_Ni,theta_0)
#'
#' @param nodei An indicator of node.
#' @param nodesLabel The vector of latent state of all genes.
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param theta_0 The value of theta_0.
#' @return The conditional probability: P(Si|S_Ni,theta_0)
#'
MRF.condProb <- function(nodei,nodesLabel,nodesNeighbor,w,W,theta_0) {
  nNodes<-length(nodesLabel)
  i.neighbor <- nodesNeighbor[[nodei]]
  h1 <- theta_0[1]
  tau11 <- theta_0[2]
  tau00 <- theta_0[3]
  if(!is.null(W)){
    i.neighbor.0 <- is.element(1:nNodes,i.neighbor)&nodesLabel==-1
    i.neighbor.1 <- is.element(1:nNodes,i.neighbor)&nodesLabel== 1
    logitp <- h1 + tau11*(sum(W[nodei,]*i.neighbor.1)) -
      tau00*(sum(W[nodei,]*i.neighbor.0))
  }else{
    i.neighbor.0 <- i.neighbor[nodesLabel[i.neighbor]==-1]
    i.neighbor.1 <- i.neighbor[nodesLabel[i.neighbor]== 1]
    logitp <- h1 + tau11*( w[nodei]*sum(nodesLabel[i.neighbor]== 1)+sum(w[ i.neighbor.1])) -
      tau00*( w[nodei]*sum(nodesLabel[i.neighbor]==-1)+sum(w[ i.neighbor.0]))
  }

  1/(1+exp(-logitp))
}

#'Calculate the summation of logP(Y|theta1)
#'
#' @param nodesLabel The vector of latent state of all genes.
#' @param theta_1 The value of theta_1.
#' @param Y The vector of DNM counts for all genes.
#' @param N The parameter in the prior for network.
#' @param mu The vector of mutability for all genes.
#' @return The summation of logP(Y|theta1)
#'
count.llk<-function(nodesLabel,theta_1,Y,N,mu){
  lambda_0<-2*N*mu
  lambda_1<-2*N*mu*exp(theta_1)
  lk<-ifelse(nodesLabel==1, dpois(Y,lambda_1),dpois(Y,lambda_0))
  sum(log(lk))
}

#'Calculate the pseudo-likelihood
#'
#' @param nodesLabel The vector of latent state of all genes.
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param theta_0 The value of theta_0.
#' @param theta_1 The value of theta_1.
#' @param Y The vector of DNM counts for all genes.
#' @param N The parameter in the prior for network.
#' @param mu The vector of mutability for all genes.
#' @return The pseudo-likelihood.
#'
pseudo.llk<-function(nodesLabel,nodesNeighbor,w,W,theta_0,theta_1,Y,N,mu)
{
  cond.lk <- unlist(lapply(1:length(nodesLabel),MRF.condProb,
                           nodesLabel=nodesLabel,nodesNeighbor=nodesNeighbor,w=w,W=W,theta_0=theta_0))
  cond.llk <-sum(log(cond.lk))
  count.llk(nodesLabel,theta_1,Y,N,mu)+cond.llk
}


#'Calculate the density for each node given its label
#'
#' @param nodesLabel The vector of latent state of all genes.
#' @param Y The vector of DNM counts for all genes.
#' @param mu The vector of mutability for all genes.
#' @param N The parameter in the prior for network.
#' @param theta_1 The value of theta_1.
#' @return The vector of density for all genes.
#'
likpw<- function(nodesLabel,Y,mu,N,theta_1=NULL) {
  lambda_0<-2*N*mu
  lambda_1<-2*N*mu*exp(theta_1)
  ifelse(nodesLabel==1,dpois(Y,lambda_1),dpois(Y,lambda_0))
}


#'Calculate the FDR for all genes
#'
#' @param Gene The vector of gene symbols.
#' @param q_0 The vector of posterior probabilities of P(S_i=-1) estimated from Gibbs.
#' @return The vector of FDR (\code{FDR}) for all genes(\code{Gene}) .
#'
Get_Post_FDR<-function(Gene,q_0){
  #Adapted from M-DATA
  #posterior should be the marginal posterior sampling from Gibbs
  Num<-length(Gene)
  lfdr<-data.frame(index=1:Num,Z_no=q_0)
  lfdr_order_no<-lfdr[order(lfdr$Z_no,decreasing = FALSE),]
  lfdr_order_no$FDR_no<-cumsum(lfdr_order_no$Z_no)/(1:Num)
  tmp_no<-lfdr_order_no[order(lfdr_order_no$index,decreasing = FALSE),]
  tmp_FDR_no<-tmp_no$FDR_no
  return(list(Gene=Gene,FDR=tmp_FDR_no))
}

#'Estimate the parameters iteratively
#'
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param Y The vector of DNM counts for all genes.
#' @param N The parameter in the prior for network.
#' @param mu The vector of mutability for all genes.
#' @param theta_0_init The initial value of theta_0.
#' @param theta_1_init The initial value of theta_1.
#' @param lambda The tuning parameter for estimation of theta_0.
#' @param threshold The threshold for convergence.
#' @param maxiter The maximum number of iteration.
#' @param initLabel The vector of initial labels for all genes.
#' @return The list containing the estimated values for \code{theta_0_est},\code{theta_1_est}
#' and the pseudo-likelihood \code{pseudollk}.
#'
icm_randsearch<-function(nodesNeighbor,w=NULL,W=NULL,
                         Y,mu,N,
                         theta_0_init,theta_1_init,
                         maxiter=100,threshold=1e-3,lambda=0.5,initLabel=NULL){
  nNodes<-length(Y)
  if(is.null(initLabel)){
    initLabel<-  2*rbinom(nNodes,1,0.5)-1
  }
  currentNodesLabel<-initLabel

  theta_0<-theta_0_init
  theta_1<-theta_1_init

  iteNum<-1;converged<-F
  #Function to update all labels in nNodes with a random sequence
  randsearch <- function(initNode,currentNodesLabel,theta_0,theta_1)
  {
    temp.label <- currentNodesLabel
    lik.1 <- likpw(rep(1,nNodes),Y,mu,N,theta_1)
    lik.0 <- likpw(rep(-1,nNodes),Y,mu,N,theta_1)
    likRatio <- lik.1/lik.0
    for (nodei in (1:nNodes)) {
      condProb<-MRF.condProb(nodei,temp.label,nodesNeighbor,w=w,W=W,theta_0=theta_0)
      temp.label[nodei]<-2*(likRatio[nodei] * (condProb /(1-condProb)) >1 )-1
      if (anyNA(temp.label)) browser();
    }
    # temp.prob<-temp.BF/(1+temp.BF)
    return(temp.label)
  }
  while (iteNum<maxiter & !converged)
  { print(paste0("ICM iteration round ",iteNum," starts!"))

    prevNodesLabel<-currentNodesLabel
    prev_theta_0<-theta_0
    prev_theta_1<-theta_1

    #algorithm step 2: update network parameter theta_0
    theta_0<-Get_theta_0(prevNodesLabel,nodesNeighbor,w,W,lambda = lambda)

    #algorithm step 3: ICM
    currentNodesLabel<-randsearch(sample(1:nNodes,1),prevNodesLabel,theta_0 = theta_0,theta_1 = prev_theta_1)

    #algorithm step 4: update data parameter theta_1
    theta_1<-Get_theta_1_Gamma(currentNodesLabel,Y,N,mu,beta_1_init = prev_theta_1)


    a<-pseudo.llk(prevNodesLabel,nodesNeighbor,w,W,prev_theta_0,prev_theta_1,Y,N,mu)
    b<-pseudo.llk(currentNodesLabel,nodesNeighbor,w,W,theta_0,theta_1,Y,N,mu)

    if(abs(a-b)/nNodes<threshold & (all(abs(theta_0-prev_theta_0)/(0.1+abs(theta_0))<threshold)) &sum(abs(theta_1-prev_theta_1))<threshold) {
      converged<-T
    }

    if(a>b & sum(abs(theta_0-prev_theta_0))<threshold) {
      print("Warning: the update does not increase conditional posterior!")
      currentNodesLabel<-prevNodesLabel
      theta_0<-prev_theta_0
      theta_1<-prev_theta_1
      converged<-T
    }
    iteNum <- iteNum+1
  }

  if(!all(converged)){print(paste0("ICM did not converge after ",iteNum," iterations!"))}else{
    print(paste0("ICM converges after ",iteNum," iterations!"))
  }
  return(list(theta_0_est=theta_0,theta_1_est=theta_1,pseudollk=b))
}


#'Calculate the Gibbs posterior
#'
#' @param theta_0 The estimated value of theta_0.
#' @param theta_1 The estimated value of theta_1.
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param Y The vector of DNM counts for all genes.
#' @param N The parameter in the prior for network.
#' @param mu The vector of mutability for all genes.
#' @param mcmc_samples The number of iterations for mcmc.
#' @param burnin The number of burn-ins for mcmc.
#' @param thin The number of thinning for mcmc.
#' @param initLabel The vector of initial labels for all genes.
#' @return The list containing \code{S_mcmc} the mcmc samples after burning-in and thinning,
#' \code{q_0} the posterior probability for calculating FDR,
#' \code{llk} the likelihood for checking convergence
#'
Gibbs.post<-function(theta_0,theta_1,nodesNeighbor,w,W=NULL,
                     Y,N,mu,
                     mcmc_samples=1010,burnin=10,thin=1,initLabel=NULL){
  nNodes <- length(Y)
  lik.1 <- likpw(rep(1,nNodes),Y,mu,N,theta_1)
  lik.0 <- likpw(rep(-1,nNodes),Y,mu,N,theta_1)
  likRatio <- lik.1/lik.0
  #Record structure
  sampling_num<-(mcmc_samples-burnin)/thin
  prob<-rep(NA,nNodes)
  llk<-rep(NA,mcmc_samples)
  S<-matrix(NA, nrow=mcmc_samples, ncol=nNodes)
  if(is.null(initLabel)){
    initLabel<-  2*rbinom(nNodes,1,0.5)-1
  }
  CurrentLabel<-initLabel
  for (i in 1:mcmc_samples) {
    # print(paste0("Sampling iteration:",i))
    if ( i %% 1000 ==0) print(paste0("mcmc iteration:",i));
    for (j in (1:nNodes)) {
      condProb<-MRF.condProb(j,CurrentLabel,nodesNeighbor,w,W,theta_0)
      prob[j]<-(likRatio[j]*condProb) / (likRatio[j]*condProb + (1-condProb))
      S[i,j]<-2*rbinom(1,1,prob[j])-1
      CurrentLabel[j]<-S[i,j]
      # browser()
    }
    llk[i]<-sum(log(prob))
  }
  S_mcmc<-S[seq(from=(burnin+1),to=mcmc_samples,by=thin),]
  q_0<-unlist(apply(S_mcmc, 2, function(x) sum(x==-1)))/sampling_num
  return(list(S_mcmc=S_mcmc,q_0=q_0,llk=llk))
}


#'Functions for sampling S_i from prior distribution
#'
#' @param true_theta_0 The true value of theta_0.
#' @param nodesNeighbor Neighbors of corresponding nodes in the network.
#' @param w The parameter in the prior for network.
#' @param W Another way to input the weight. For internal use.
#' @param mcmc_samples The number of iterations for mcmc.
#' @param burnin The number of burn-ins for mcmc.
#' @param S_init The vector of initial labels for all genes.
#' @return The list containing \code{S} all mcmc samples,
#' \code{S_mcmc} the mcmc samples after burning-in,
#' \code{llk} the likelihood for checking convergence
#'
Gibbs.trueParam<-function(true_theta_0,nodesNeighbor,w,W=NULL,mcmc_samples=1000,burnin=500,S_init=NULL){
  nNodes <- length(w)
  #Record structure
  condProb<-rep(NA,nNodes)
  llk<-rep(NA,mcmc_samples)
  S<-matrix(NA, nrow=mcmc_samples, ncol=nNodes)
  CurrentLabel<-2*rbinom(nNodes,1,0.5)-1
  for(i in 1:mcmc_samples){
    for(j in (1:nNodes)){
      condProb[j]<-MRF.condProb(j,CurrentLabel,nodesNeighbor,w,W,true_theta_0)
      S[i,j]<-2*rbinom(1,1,condProb[j])-1
      CurrentLabel[j]<-S[i,j]
    }
    llk[i]<-sum(log(condProb))
  }
  # S_mcmc<-unlist(apply(S, 1, getmode))
  S_mcmc<-S[(burnin+1):mcmc_samples,]
  return(list(S=S,S_mcmc=S_mcmc,llk=llk))
}

icm_randsearch_theta0<-function(nodesNeighbor,w=NULL,W=NULL,
                                Y,mu,N,
                                theta_0_init,theta_1_init,
                                maxiter=100,threshold=1e-3,lambda=0.5,initLabel=NULL){
  nNodes<-length(Y)
  if(is.null(initLabel)){
    initLabel<-  2*rbinom(nNodes,1,0.5)-1
  }
  currentNodesLabel<-initLabel

  theta_0<-theta_0_init
  theta_1<-theta_1_init

  iteNum<-1;converged<-F
  #Function to update all labels in nNodes with a random sequence
  randsearch <- function(initNode,currentNodesLabel,theta_0,theta_1)
  {
    temp.label <- currentNodesLabel
    lik.1 <- likpw(rep(1,nNodes),Y,mu,N,theta_1)
    lik.0 <- likpw(rep(-1,nNodes),Y,mu,N,theta_1)
    likRatio <- lik.1/lik.0
    for (nodei in (1:nNodes)) {
      condProb<-MRF.condProb(nodei,temp.label,nodesNeighbor,w=w,W=W,theta_0=theta_0)
      temp.label[nodei]<-2*(likRatio[nodei] * (condProb /(1-condProb)) >1 )-1
      if (anyNA(temp.label)) browser();
    }
    # temp.prob<-temp.BF/(1+temp.BF)
    return(temp.label)
  }
  while (iteNum<maxiter & !converged)
  { print(paste0("ICM iteration round ",iteNum," starts!"))

    prevNodesLabel<-currentNodesLabel
    prev_theta_0<-theta_0
    # prev_theta_1<-theta_1

    #algorithm step 2: update network parameter theta_0
    theta_0<-Get_theta_0(prevNodesLabel,nodesNeighbor,w,W,lambda = lambda)

    #algorithm step 3: ICM
    currentNodesLabel<-randsearch(sample(1:nNodes,1),prevNodesLabel,theta_0 = theta_0,theta_1 =theta_1)

    a<-pseudo.llk(prevNodesLabel,nodesNeighbor,w,W,prev_theta_0,theta_1,Y,N,mu)
    b<-pseudo.llk(currentNodesLabel,nodesNeighbor,w,W,theta_0,theta_1,Y,N,mu)

    if(abs(a-b)/nNodes<threshold & (all(abs(theta_0-prev_theta_0)/(0.1+abs(theta_0))<threshold))<threshold) {
      converged<-T
    }

    if(a>b & sum(abs(theta_0-prev_theta_0))<threshold) {
      print("Warning: the update does not increase conditional posterior!")
      currentNodesLabel<-prevNodesLabel
      theta_0<-prev_theta_0
      # theta_1<-prev_theta_1
      converged<-T
    }
    iteNum <- iteNum+1
  }

  if(!all(converged)){print(paste0("ICM did not converge after ",iteNum," iterations!"))}else{
    print(paste0("ICM converges after ",iteNum," iterations!"))
  }
  return(list(theta_0_est=theta_0,theta_1_est=theta_1,pseudollk=b,icm=currentNodesLabel))
}
