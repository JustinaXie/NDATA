## Main function for HMRF analysis: perform HMRF analysis to find risk genes from the network estimated by PNS algorithm
## Input:
##      pv - vector, p-values of genes in the network
##      graph - matrix, estimated network in matrix format
##      covGene - vector, names of genes in the additional covariate of HMRF analysis
##      fixedGene - vector, names of genes to have fixed hidden states 1 in HMRF analysis
##      pthres_hmrf - threshold for p-values in HMRF initialization, default=0.05
##      iter - number of iteration, default=100
##      trimthres - threshold for Z score trimming, default=5
## Output:
##      result_DAWN - a list object representing the estimated statistics of genes
##          result_DAWN$Iupdate - updated I
##          result_DAWN$post - posterior distribution of I
##          result_DAWN$b0 - estimated b
##          result_DAWN$b1 - estimated c
##          result_DAWN$b2 - estimated coefficient for additional covariate (default=0 if additional covariate is not included)
##          result_DAWN$mu1 - estimated mean
##          result_DAWN$sigmas - estimated variance
DAWN_HMRF <- function(pv, graph,covGene, fixedGene, pthres_hmrf=0.05, iter=100, trimthres=5,initLabel) {
  ## initialization
  zscore <- qnorm(1 - pv)
  if (missing(pthres_hmrf)) {
    pthres_hmrf <- 0.05
  }
  if (missing(iter)) {
    iter <- 100
  }
  if (missing(trimthres)) {
    trimthres <- 4
  }
  addCovGene <- hasArg(covindex)
  addFixedGene <- hasArg(fixindex)
  if (addCovGene) {
    covindex <- rep(0, length(genelist_PNS))
    covindex[(genelist_PNS) %in% covGene] <- 1
  }
  if (addFixedGene) {
    fixindex <- rep(0, length(genelist_PNS))
    fixindex[(genelist_PNS) %in% fixedGene] <- 1
    fixindex <- (zscore > trimthres) | fixindex + 0
  } else {
    fixindex <- (zscore > trimthres) + 0
  }

  Istart<-(pv < pthres_hmrf) + 0
  Istart[which(fixindex==1)] <- 1
  b1_1 <- 0
  b0_1 <- 0
  b2_1 <- 0
  Iupdate <- Istart
  mu1 <- mean(zscore[Istart == 1 & fixindex == 0])
  sigmas1 <- (sd(zscore[Istart == 0]))^2
  sigmas2 <- (sd(zscore[Istart == 1 & fixindex == 0]))^2
  posterior <- rep(0, length(zscore))

  for (iteri in 1:iter) {
    ## Algorithm2, step2(a) optimize b, c (and d)
    print(paste0("iteration: ",iteri," begins!"))
    res <- optimize_params(graph, Iupdate, addCovGene, covindex, 20)
    b0 <- res[1]
    b1 <- res[2]
    b2 <- res[3]
    if (addCovGene) {
      cat("iteration", iteri-1, ": b0 =", b0_1,", b1 =", b1_1, ", b2 =", b2_1, ", mu1 =", mu1, ", sigmas1 =", sigmas1, "\n")
    } else {
      cat("iteration", iteri-1, ": b0 =", b0_1, ", b1 =", b1_1, ", mu1 =", mu1, ", sigmas1 =", sigmas1, "\n")
    }
    if (abs(b1_1 - b1) < 0.001 & abs(b0_1 - b0) < 0.001 & abs(b2_1 - b2) < 1e-05) {
      break
    }
    b1_1 <- b1
    b0_1 <- b0
    b2_1 <- b2
    ## Algorithm2, step2(b) update I
    Ibefore = Iupdate
    for (i in 1:length(pv)) {
      new_vec <- compute_new(i, b0, b1, b2, Iupdate, graph[i, ], addCovGene, covindex[i])
      new1 <- new_vec[1]
      new2 <- new_vec[2]
      p1 <- dnorm(zscore[i], mu1 * Iupdate[i], sqrt(sigmas2 * Iupdate[i] + sigmas1 * (1 - Iupdate[i]))) /(1 +
                                                                                                            exp(new2-new1))  ##sigmas2 is the alternative, sigmas1 is the null var
      p2 <- dnorm(zscore[i], mu1 * (1 - Iupdate[i]), sqrt(sigmas2 * (1 - Iupdate[i]) + sigmas1 * Iupdate[i])) /
        (1 + exp(new1-new2))
      if(is.na(p1)|is.na(p2)) browser();

      if (Iupdate[i] == 1)
        posterior[i] <- p1/(p1 + p2)
      if (Iupdate[i] == 0)
        posterior[i] <- p2/(p1 + p2)
      if (p2 > p1) {
        Iupdate[i] <- 1 - Iupdate[i]
      }
      if (fixindex[i] != 0) {
        Iupdate[i] <- 1
      }
    }
    ## Algorithm2, step2(c) update parameters
    mu1 <- sum(posterior[fixindex == 0] * zscore[fixindex == 0])/sum(posterior[fixindex == 0])
    sigmas2 <- sum(posterior[fixindex == 0] * (zscore[fixindex == 0] - mu1)^2)/sum(posterior[fixindex == 0])
    sigmas1 <- sum((1 - posterior[fixindex == 0]) * (zscore[fixindex == 0])^2)/sum(1 - posterior[fixindex == 0])
    sigmas <- (sigmas1 * sum(posterior[fixindex == 0]) + sigmas2 * sum(1 - posterior[fixindex == 0]))/length(posterior)
    sigmas2 <- sigmas
    sigmas1 <- sigmas
  }
  result_DAWN <- list()
  result_DAWN$Iupdate <- Iupdate
  result_DAWN$post <- posterior
  result_DAWN$b0 <- b0
  result_DAWN$b1 <- b1
  result_DAWN$b2 <- b2
  result_DAWN$mu1 <- mu1
  result_DAWN$sigmas <- sigmas

  if (!check_b0(I = Iupdate, b0 = b0)) {
    cat('\n\nWARNING: DAWN identified a large number of risk genes.
          Assumptions of the model may be false.
          The set of risk genes likely contains many false positives.\n')
  }
  if(is.na(b0)) browser();
  if(is.na(b1)) browser();
  if(anyNA(Iupdate)) browser();
  if(is.na(check_b1(G = graph, I = Iupdate, b0 = b0, b1 = b1))) print(paste0("Debug:",Iupdate%*%graph%*%Iupdate))
  if(is.na(check_b1(G = graph, I = Iupdate, b0 = b0, b1 = b1))) print(paste0("Debug2:",sum(Iupdate)))
  if (!check_b1(G = graph, I = Iupdate, b0 = b0, b1 = b1)) {
    cat('\n\nWARNING: Weak connectivity among risk genes in the input graph.
          Assumptions of the model appear to be false.
          The set of risk genes likely contains many false positives.\n')
  }
  return(result_DAWN)
}


##  Formating final results
##  Input:
##      finalposter - final posterior of hidden states I
##      kgene - genes in the graph
##      pv - p-values of genes in the graph
## Output:
##      report - a data frame containing information for genes and their corresponding posteriors, GLA p-values, and FDR
result_report <- function(finalposter, genes, pv) {
  usepost <- finalposter
  rankpost <- sort(usepost)
  localfdr <- rep(0, length(usepost))

  for (i in 1:length(localfdr)) {
    localfdr[i] <- mean(rankpost[1:i])
  }

  flocalfdr <- rep(0, length(localfdr))
  rankp <- rank(usepost, ties.method = "random")
  flocalfdr <- localfdr[rankp]

  report <- data.frame(genes, usepost, pv, flocalfdr)
  names(report) <- c("gene", "posterior", "GLA.pvalue", "FDR")
  return(report)
}


####### supporting functions for HMRF analysis #######

## Compute intermediate results for Algorithm2, step2(a)
## Input:
##      i - index of current gene
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      b2 - parameter d in Ising model (coefficient for additional covariate, default=0 if additional covariate is not used)
##      Iupdate - hidden states I
##      graph_i - connectivity of current gene to other genes in PNS (i'th row in PNS graph)
##      addCovGene - whether to use covGene information or not
##      covindex_i - if current gene is in covindex
## Output:
##      a vector of intermediate results (new1, new2)
compute_new <- function(i, b0, b1, b2, Iupdate, graph_i, addCovGene, covindex_i) {
  if (addCovGene) {
    new1 <- (b0 * Iupdate[i] + b1 * Iupdate[i] * t(graph_i) %*% Iupdate + b2 * Iupdate[i] * covindex_i)
    new2 <- (b0 * (1 - Iupdate[i]) + b1 * (1 - Iupdate[i]) * t(graph_i) %*% Iupdate + b2 * (1 - Iupdate[i]) *
               covindex_i)
  } else {
    new1 <- (b0 * Iupdate[i] + b1 * Iupdate[i] * t(graph_i) %*% Iupdate)
    new2 <- (b0 * (1 - Iupdate[i]) + b1 * (1 - Iupdate[i]) * t(graph_i) %*% Iupdate)
  }
  return(c(new1, new2))
}

## Estimate parameters in Ising model by maximizing the pseudo likelihood
## Input:
##      graph - graph estimated by PNS algorithm
##      Iupdate - hidden indicator vector I
##      addCovGene - whether to use covGene information or not
##      covindex - vector, indicates whether a gene is in the addtional covariate
##      times - number of iterations
## Output:
##      a vector of estimated paramter values - c(b0, b1, b2)
##        b0 - parameter b in Ising model
##        b1 - parameter c in Ising model
##        b2 - parameter d in Ising model (coefficient for additional covariate, default=0 if additional covariate is not used)
optimize_params <- function(graph1, Iupdate, addCovGene, covindex, times) {
  ## initialization
  b0_iter <- 0
  b1_iter <- 0
  b2_iter <- 0
  risknode.num <- Iupdate %*% graph1
  if (addCovGene) {
    for (k in 1:times) {
      b0_iter_1 <- optimize(plf_covGene, c(-20, 0), b1 = b1_iter, b2 = b2_iter, risknode.num = risknode.num, zupdate = Iupdate,
                            covindex = covindex, maximum = T)$maximum
      b1_iter_1 <- optimize(plf_covGene, c(0, 10), b0 = b0_iter_1, b2 = b2_iter, risknode.num = risknode.num, zupdate = Iupdate,
                            covindex = covindex, maximum = T)$maximum
      b2_iter_1 <- optimize(plf_covGene, c(0, 10), b1 = b1_iter_1, b0 = b0_iter_1, risknode.num = risknode.num, zupdate = Iupdate,
                            covindex = covindex, maximum = T)$maximum
      if (abs(b0_iter_1 - b0_iter) < 10^-5 & abs(b1_iter_1 - b1_iter) < 10^-5 & abs(b2_iter_1 - b2_iter) < 10^-5) {
        break
      }
      b0_iter <- b0_iter_1
      b1_iter <- b1_iter_1
      b2_iter <- b2_iter_1
    }
  } else {
    for (k in 1:times) {
      b0_iter_1 <- optimize(plf, c(-20, 0), b1 = b1_iter, risknode.num = risknode.num, zupdate = Iupdate, maximum = T)$maximum
      b1_iter_1 <- optimize(plf, c(0, 10), b0 = b0_iter_1, risknode.num = risknode.num, zupdate = Iupdate, maximum = T)$maximum
      if (abs(b0_iter_1 - b0_iter) < 10^-5 & abs(b1_iter_1 - b1_iter) < 10^-5) {
        break
      }
      b0_iter <- b0_iter_1
      b1_iter <- b1_iter_1
    }
  }

  return(c(b0_iter, b1_iter, b2_iter))
}

## Compute pseudo likelihood
## Input:
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      risknode.num - number of risk nodes connected to each node in the graph estimated by PNS algorithm
##      zupdate - hidden indicators I
## Output:
##      fvalue - pseudo likelihood value
plf <- function(b0, b1, risknode.num, zupdate) {
  fv.fun <- function(i) {
    new1 <- exp(b0 * zupdate[i] + b1 * zupdate[i] * risknode.num[i])
    new2 <- exp(b0 * (1 - zupdate[i]) + b1 * (1 - zupdate[i]) * risknode.num[i])
    fvalue <- log(new1/(new1 + new2))
  }
  fvalue <- sum(sapply(c(1:length(zupdate)), fv.fun))
  return(fvalue)
}

## Compute pseudo likelihood with additional covariate added
## Input:
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      b2 - parameter d in Ising model (coefficient for covGenes)
##      risknode.num - number of risk nodes connected to each node in the graph estimated by PNS algorithm
##      covindex - vector, indicates whether a gene is in the addtional covariate
##      zupdate - hidden indicator vector I
## Output:
##      fvalue - pseudo likelihood value
plf_covGene <- function(b0, b1, b2, risknode.num, covindex, zupdate) {
  fv.fun <- function(i) {
    new1 <- exp(b0 * zupdate[i] + b1 * zupdate[i] * risknode.num[i] + b2 * covindex[i] * zupdate[i])
    new2 <- exp(b0 * (1 - zupdate[i]) + b1 * (1 - zupdate[i]) * risknode.num[i] + b2 * covindex[i] * (1 - zupdate[i]))
    fvalue <- log(new1/(new1 + new2))
  }
  fvalue <- sum(sapply(c(1:length(zupdate)), fv.fun))
  return(fvalue)
}

## Check if b0 value is reasonable
## Input:
##      I - estimated hidden indicator vector
##      b0 - estimated b0
##      thres_b0 - threshold for evaluating b0
## Output:
##      pass - a boolean variable, True if b0 passes the test, False otherwise
check_b0 <- function(I, b0, thres_b0 = 0.05){
  ## Check if P(sum(I)>thres) < thres under the distribution that:
  ## P(sum(I)) ~ binom(n,exp(b0)/(1+exp(b0)))
  n <- length(I)
  q <- floor(n * thres_b0);
  p <- exp(b0) / (1 + exp(b0))
  tail_prob <- pbinom(q = q, size=n, prob = p, lower.tail = F) ##P(sum(I)>thres)
  pass <- (tail_prob < thres_b0)
  return(pass)
}

## Check if b1 value is reasonable
## Input:
##      G - adjacency matrix for graph
##      I - estimated hidden indicator vector
##      b0 - estimated b0
##      b1 - estimated b1
##      thres_b1 - threshold for evaluating b1
## Output:
##      pass - a boolean variable, True if b0 passes the test, False otherwise
check_b1 <- function(G, I, b0, b1, thres_b1 = 0.1) {
  ## Check if b1*I'GI is significant comparing to b0*I
  pass <- abs(b1*I%*%G%*%I / (b0*sum(I)))[1,1] > thres_b1
  return(pass)
}


