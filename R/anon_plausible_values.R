#############################################################
# Functions relating to plausible values
#############################################################

## Plausible Values
# pv1 <- function(b,a,first,last,score,mu=0,sigma=2)
# {
#   tmp=.C("PV1",
#          as.double(b),theta_
#          as.integer(a),
#          as.integer(first-1),
#          as.integer(last-1),
#          as.double(mu),
#          as.double(sigma),
#          as.integer(score),
#          as.integer(rep(0,length(score))),
#          as.integer(length(score)),
#          as.integer(length(first)),
#          as.integer(1),
#          as.double(0*score))[[12]]
#   return(tmp)
# }

# Wrapper for C function that generates plausible values using a straightforward rejection algorithm.
# This one does not use recycling
# Normal prior(s)
# @param b          vector of beta's per item_score, including 0 category, ordered by item_id, item_score
# @param a          vector of discriminations per item_score, inclusing 0 category, ordered by item_id, item_score
# @param score      vector of sum scores
# @param npv        nbr of plausible values per sum score
# @param pop        vector of indices of population for each sum score
# @param mu, sigma  prior mu and sigma for each population
# @return
# matrix with in each row npv plausible values for each sum score
pv_ <- function(b,a,first,last,score,npv,mu,sigma,pop)
{
  tmp=.C("PV",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(mu),
         as.double(sigma),
         as.integer(score),
         as.integer(pop-1),
         as.integer(length(score)),
         as.integer(length(first)),
         as.integer(npv),
         as.double(rep(0*score,npv)))[[12]]
  tmp=as.vector(tmp)
  if (npv>1) dim(tmp)=c(length(score),npv)
  return(tmp)
}

# score is a vector of scores
# mu and sigma are scalar. One prior assumed
# returned is a matrix with for each score (=row) npv plausible values
pv_recycle <- function(b,a,first,last,score,npv,mu,sigma)
{
  # could be made slightly faster if score already sorted
  #stopifnot(npv==1)
  #pop = unique(pop)
  #stopifnot(length(pop)==1)
  
  if (sigma<0) stop('prior standard deviation must be positive')
  mx = sum(a[last])
  
  scrs = tibble(score=score) %>% 
    mutate(indx = row_number()) %>%
    arrange(.data$score)
  
  scr_tab = scrs %>%
    group_by(.data$score) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    right_join(tibble(score=0L:mx),by='score') %>%
    mutate(n=coalesce(n,0L)) %>%
    arrange(.data$score) %>%
    pull(.data$n) * npv
  
  cscr_tab = cumsum(scr_tab)
  
  theta = rep(0,length(score)*npv)
  
  
  tmp=.C("PVrecycle",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(mu),
         as.double(sigma),
         as.integer(scr_tab),
         as.integer(cscr_tab),
         as.integer(length(score)*npv),
         as.integer(length(first)),
         as.integer(max(a)),
         as.double(theta))[[12]]
  
  if(npv>1)
  {
    tmp = matrix(tmp, length(score), npv,byrow=TRUE)[order(scrs$indx),,drop=FALSE]
  } else
  {
    tmp = tmp[order(scrs$indx)]
  }
  
  return(tmp)
}


# Wrapper for C function that generates plausible values using a straightforward rejection algorithm
# Mixture of two normals as flexible prior
# @param b          vector of beta's per item_score, including 0 category, ordered by item_id, item_score
# @param a          vector of discriminations per item_score, inclusing 0 category, ordered by item_id, item_score
# @param score      vector of sum scores
# @param npv        nbr of plausible values per sum score
# @param p          prior group membership probabilities
# @param mu, sigma  mu and sigma for each population
# @param pop        vector of indices of population for each sum score. Currently not used
# @return
# matrix with in each row npv plausible values for each sum score
pv_mix <- function(b,a,first,last,score,npv,p,mu,sigma, pop)
{
  tmp=.C("PVMix",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(p),
         as.double(mu),
         as.double(sigma),
         as.integer(score),
         as.integer(length(score)),
         as.integer(length(first)),
         as.integer(npv),
         as.double(rep(0*score,npv)))[[12]]
  tmp=as.vector(tmp)
  if (npv>1) dim(tmp)=c(length(score),npv)
  return(tmp)
}

## Uses recycling to produce npv plausible values per score in scores
#  on a test defined by first and last. Scores defined to be 
#  0 to the maxScore on the test
#  If a prior sample is provided this is used; Otherwise prior is normal(mu,sigma)
#  A is a vector like a with alternative weightes. Used in theta_tables.
#  nscore is optional: a vector containing for each possible test score a frequency
#  is !is.null(nscore) and !is.null(A) I assume that nscore refers to scores with A as weights
#TO DO: This function needs to go in favour of pv_recycle
recycle_pv = function(b, a, first, last, npv=1, mu=0, sigma=2, nscore = NULL, prior_sample=NULL, A=NULL)
{
  if (is.null(A)){
    scores=0:sum(a[last])
    not_possible=setdiff(scores,possible_scores(a,first,last))
  }else {
    if (length(A)!=length(a)) stop("Wrong vector A in recycle_pv")
    if (any(A)<0) stop("Negative weights in vector A in recycle_pv")
    scores=0:sum(A[last])
    not_possible=setdiff(scores,possible_scores(A,first,last))
  }
  
  if (is.null(nscore)) nscore=rep(1,length(scores))
  n=nscore*npv
  if (length(not_possible)>0) n[not_possible+1]=0
  R=rep(0,sum(n))
  
  nI=length(first)
  if (is.null(A))
  {
    if (is.null(prior_sample))
    {
      if (sigma<0) stop('prior standard deviation must be positive')
      tmp=.C("recyclePV",
             as.double(b),
             as.integer(a),
             as.integer(first-1),
             as.integer(last-1),
             as.integer(nI),
             as.integer(n),
             as.double(mu),
             as.double(sigma),
             as.double(R))[[9]]
    }else
    {
      nprior=length(prior_sample)
      if (nprior<1) stop("wrong prior sample")
      tmp=.C("recyclePV2",
             as.double(b),
             as.integer(a),
             as.integer(first-1),
             as.integer(last-1),
             as.integer(nI),
             as.integer(n),
             as.double(prior_sample),
             as.double(nprior),
             as.double(R))[[9]]
    }
  }else
  {
    if (sigma<0) stop('prior standard deviation must be positive')
    tmp=.C("recyclePVaA",
           as.double(b),
           as.integer(a),
           as.integer(A),
           as.integer(first-1),
           as.integer(last-1),
           as.integer(nI),
           as.integer(n),
           as.double(mu),
           as.double(sigma),
           as.double(R))[[10]]
  }
  
  
  dim(tmp)=c(npv,length(scores))
  tmp=t(tmp)
  tmp[not_possible+1,]=NA
  
  return(tmp)
}

####################### FUNCTIONS TO UPDATE PRIOR of PLAUSIBLE VALUES

# Given samples of plausible values from one or more normal distributions
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values
#
# @param pv     vector of plausible values
# @param pop    vector of population indexes. These indices refer to mu and sigma
# @param mu     current means of each population
# @param sigma  current standard deviation of each population
#
# @return       a sample of mu and sigma

#TO DO: Adapt when there are multiple groups. Use hyperprior to keep means together
update_pv_prior<-function(pv, pop, mu, sigma)
{
  min_n=5 # no need to update variance when groups are very small
  for (j in 1:length(unique(pop)))
  {
    pv_group=subset(pv,pop==j)
    m=length(pv_group)
    if (m>min_n)
    {
      pv_var=var(pv_group)         
      sigma[j] = sqrt(1/rgamma(1,shape=(m-1)/2,rate=((m-1)/2)*pv_var))
    }
    pv_mean=mean(pv_group)
    mu[j] = rnorm(1,pv_mean,sigma[j]/sqrt(m)) 
  }
  return(list(mu=mu, sigma=sigma))
}


# Given samples of plausible values from a mixture of two normal distributions
# this function samples means and variances of plausible values from their posterior
# It updates the prior used in sampling plausible values. 
#
# Loosely based on Chapter 6 of 
# Marin, J-M and Robert, Ch, P. (2014). Bayesian essentials with R. 2nd Edition. Springer: New-York
# but works surpringly better than the official implementation in bayess
#
# @param pv     vector of plausible values
# @param p  vector of group membership probabilities
# @param mu     current means of each group
# @param sigma  current standard deviation of each group
#
# @return       a sample of p, mu and sigma
update_pv_prior_mixnorm = function (pv, p, mu, sigma) {
  n = length(pv)
  z = rep(0, n)
  nj = c(0,0); 
  l = rep(1,2)
  v = rep(5,2)
  # Burnin and nr of Iterations with 0<Bin<nIter
  Bin = 1 
  nIter = 2
  
  mug = sig = pgr = matrix(0,nIter,2)
  mug[1,] = mu
  sig[1,] = sigma
  pgr[1,] = p
  mean_pv = mean(pv)
  var_pv = var(pv)
  for (i in 1:nIter)
  {
    ## latent group membership
    for (t in 1:n)
    {
      prob = pgr[i,]*dnorm(pv[t], mean = mug[i,], sd = sig[i,])
      if (sum(prob)==0) prob=c(0.5,0.5)
      z[t] = sample(1:2, size = 1, prob=prob)
    }
    ## Means
    for (j in 1:2)
    {
      nj[j]  = sum(z==j) 
      mug[i,j]  = rnorm(1, mean = (l[j]*mug[i,j]+nj[j]*mean(pv[z==j]))/(l[j]+nj[j]), 
                        sd = sqrt(sig[i,j]^2/(nj[j]+l[j])))
    }
    ## Vars 4=var_pv
    for (j in 1:2) 
    {
      if (nj[j]>1)
      {
        var_pvj = var(pv[z==j])
        sig[i,j] = sqrt(1/rgamma(1, shape = 0.5*(v[j]+nj[j]) ,
                                 rate = 0.5*(var_pv + nj[j]*var_pvj + 
                                               (l[j]*nj[j]/(l[j]+nj[j]))*(mean_pv - mean(pv[z==j]))^2)))
      }
    }
    ## membership probabilities
    p = rgamma(2, shape = nj + 1, scale = 1)
    pgr[i,] = p/sum(p)
  }
  return(list(p = colMeans(pgr[Bin:nIter,]), mu = colMeans(mug[Bin:nIter,]), sigma = colMeans(sig[Bin:nIter,]), grp = z))
}


# Plausible values. This one uses recycling 
# @param x                tibble(booklet_id <char or int>, person_id <char>, sumScore <int>, pop <int>)
# @param design           list: names: as.character(booklet_id), values: tibble(first <int>, last <int>) ordered by first
# @param b                vector of b's per item_score, including 0 category, ordered by item_id, item_score
# @param a                vector of weights per item_score, inclusing 0 category, ordered by item_id, item_score
# @param nPV              number of plausible values to return per person
### If the parameters are estimated with calibrate_Bayes:
# @param from             burn-in: first sample of b that is used
# @param by               every by-th sample of b is used. If by>1 auto-correlation is avoided.
# @param prior.dist       Prior distribution
#
# @return
# tibble(booklet_id <char or int>, person_id <char>, sumScore <int>, nPV nameless columns with plausible values)
pv = function(x, design, b, a, nPV, from = 20, by = 5, prior.dist = c("normal", "mixture"))
{
  prior.dist = match.arg(prior.dist)
  nPop = length(unique(x$pop))
  n_prior_updates = 10
  
  x$booklet_id = as.character(x$booklet_id)
  
  #if (prior.dist == "mixture") priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2), grp=sample(1:2, length(x$pop), replace=T, prob=c(0.5,0.5)))
  if (prior.dist == "mixture")
  {
    priors = list(p=c(0.6,0.4), mu=c(0,0.1), sigma=c(2,2))
    x$grp = sample(1:2, nrow(x), replace=T, prob=c(0.5,0.5))
  } else if (prior.dist == "normal")
  {
    priors = list(mu=rep(0,nPop), sigma=rep(4,nPop))
  }
  
  
  if (is.matrix(b))
  {
    which.pv = seq(from,(from-by)*(from>by)+by*nPV,by=by)
    nIter=max(which.pv)
    if (nrow(b)<nIter){
      stop(paste("at least", as.character(nIter), "samples of item parameters needed in function pv"))
    }
    out_pv=matrix(0,length(x$sumScore),nPV)
    
    apv=1
    pb = txtProgressBar(min=0, max=nIter)
    for(iter in 1:nIter)
    {
      if (prior.dist == "mixture")
      {
        x = x %>% 
          group_by(.data$booklet_id, .data$grp) %>%
          mutate(PVX = pv_recycle(b[iter,], a, 
                                  design[[.data$booklet_id[1]]]$first, 
                                  design[[.data$booklet_id[1]]]$last, 
                                  .data$sumScore, 1, 
                                  priors$mu[.data$grp[1]], priors$sigma[.data$grp[1]])) %>%
          ungroup() 
        
        priors = update_pv_prior_mixnorm(x$PVX, priors$p, priors$mu, priors$sigma)
        x$grp = priors$grp
        
      } else if (prior.dist == "normal") 
      {
        x = x %>% 
          group_by(.data$booklet_id,.data$pop) %>%
          mutate(PVX = pv_recycle(b[iter,], a, 
                                  design[[.data$booklet_id[1]]]$first, 
                                  design[[.data$booklet_id[1]]]$last, 
                                  .data$sumScore, 1, 
                                  priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
          ungroup()
        
        priors = update_pv_prior(x$PVX, x$pop, priors$mu, priors$sigma)
      }
      
      if (iter == which.pv[apv])
      {
        colnames(x)[colnames(x)=='PVX'] = paste0('PV', iter)
        apv=apv+1
      }
      setTxtProgressBar(pb, value=iter)
    }
    close(pb)
    return( select(x, .data$booklet_id, .data$person_id,.data$sumScore, matches('PV\\d+')))
    
  }else
  {
    for(iter in 1:n_prior_updates) 
    {
      
      pv = x %>%
        group_by(.data$booklet_id,.data$pop) %>%
        mutate(pv = pv_recycle(b, a, 
                               design[[.data$booklet_id[1]]]$first, 
                               design[[.data$booklet_id[1]]]$last, 
                               .data$sumScore, 1, 
                               priors$mu[.data$pop[1]], priors$sigma[.data$pop[1]])) %>%
        ungroup() 
      
      priors = update_pv_prior(pv$pv,pv$pop,priors$mu,priors$sigma)
    }
    
    return(  
      x %>% 
        group_by(.data$booklet_id, .data$pop) %>%
        do({
          bkID = .$booklet_id[1]
          popnbr = .data$pop[1]
          out_pv = pv_recycle(b, a, design[[bkID]]$first, design[[bkID]]$last, .$sumScore, 
                              nPV, priors$mu[popnbr], priors$sigma[popnbr])
          data.frame(.$person_id, .$sumScore, as.data.frame(out_pv), stringsAsFactors = FALSE)
        }) %>%
        ungroup() %>%
        select(-.data$pop)) 
  }
}
