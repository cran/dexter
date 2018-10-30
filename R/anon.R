
# simulate 0-based index of responses to a single item. Adapted for inclusion zero
# the index i is local in that it refers to an item-index in first and last
renorm <- function(b,a,theta,first,last,i)
{
  m=length(theta)
  x=rep(0,m)
  if ((i<0)|(i>length(first))) stop("wrong item-indx in renorm")
  tmp=.C("sampleNRM",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(i-1),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(m),
         as.integer(x))[[8]]
  return(tmp)
}

# simulate responses to a single item. 
# NOt adapted for inclusion of parameter for 
# zero category
renorm0 <- function(b,a,theta,first,last,i)
{
  m=length(theta)
  x=rep(0,m)
  tmp=.C("sampleNRM0",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(i-1),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(m),
         as.integer(x))[[8]]
  return(tmp)
}

#void sampleNRM2(double *theta, double *b, int *a, int *first, int *last, int *nI, int *m, int *test_score)
# simulate test-scores rather then response patterns. Adapted for inclusion zero
rscore <- function(theta,b,a,first,last, cntr=NULL, use_b_matrix=FALSE)
{
  if(use_b_matrix)
  {
    if(is.null(cntr)) stop('use_b_matrix is true, need a counter')
    b = b[cntr$get(),]
  }
  m=length(theta)
  nI=length(first)
  x=rep(0,m)
  mxa = max(a[last])
  tmp=.C("sampleNRM2",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nI),
         as.integer(m),
         as.integer(mxa),
         as.integer(x))[[9]]
  return(tmp)
}


# Expected scores given a single ability value theta
E_score=function(theta,b,a,first,last)
{
  escore=as.double(0.0)
  n=length(first)
  tmp=.C("Escore",
         as.double(theta),
         as.double(escore),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(n))[[2]]
  return(tmp)
}



# Expected distribution given one ability theta
pscore <- function(theta, b, a, first, last)
{
  g=elsym(b, a, first, last)
  p=rep(0, length(g))
  for (s in 1:length(g))
  {
    p[s]=g[s]*exp((s-1)*theta)
  }
  return(p/sum(p))
}

####################################################
# Computes likelihood and test information for internal use
#
# For a vector of thetas it returns:
# l = a matrix (n of response cats * length of theta) of the likelihood or log-likelihood if log=TRUE
# I = a vector of the information function computed at each theta = sum(P'^2/PQ)
# J = something sum(P'P"/PQ) 
# The vector theta can be a set of quadrature points or 
# estimates to compute their SE
#
# Note: can not deal with Inf or NA values in theta
IJ_ = function(b,a,first, last, theta, log=FALSE)
{
  nI = length(first)
  nT = length(theta)
  max_ncat = max(last-first) + 1L
  I = rep(0, nI * nT)
  J = rep(0, nI * nT)
  logFi = rep(0, nT)
  
  tmp=.C("IJ_c",
         as.double(theta),
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.double(I),
         as.double(J),
         as.double(logFi),
         as.integer(nI),
         as.integer(nT),
         as.integer(max_ncat))

  scores = 0:sum(a[last])
  l = sweep(outer(scores,theta), 2, tmp[[8]], '-')
  if (!log) l=exp(l)
  return(list(I=colSums(matrix(tmp[[6]],nI,nT)), J=colSums(matrix(tmp[[7]],nI,nT)), l=l))
}


#### ML estimation of theta
# Now uses false position to determine the MLE
# se is the option to calculate standard-errors or not
theta_MLE <- function(b,a,first,last, se=FALSE)
{
  mx=sum(a[last])
  theta=rep(0,mx-1)
  max_a = max(a)

  n=length(first)
  theta = .C("theta_mle_fp",
           as.double(theta),
           as.double(b),
           as.integer(a),
           as.integer(first-1L),
           as.integer(last-1L),
           as.integer(n),
           as.integer(mx),
           as.integer(max_a))[[1]]

  sem = NULL
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem = c(NA, 1/sqrt(f$I), NA)
  }
  theta=c(-Inf,theta,Inf)
  
  return(list(theta=theta,se=sem))
}



##### EAP based on npv (default 500) plausible values ####
# Uses recycling to get npv plausible values for each sum score.
# Options:
# Smoothing is an option to avoid non-monotonicity due to random fluctuations
# se is the option to calculate standard-errors or not
### In R-code:
# n=rep(npv,length(score))
# R=matrix(0,length(score),npv)
# while (any(n>0))
# {
#   atheta=rnorm(1,mu,sigma)
#   sm=rscore(b,a,atheta,first,last)+1
#   if (n[sm]>0)
#   {
#     R[sm,n[sm]]=atheta
#     n[sm]=n[sm]-1
#   }
# }
# The final argument allows EAPs based on A-weighted score, where A need not equal a.
####
theta_EAP <- function(b, a, first, last, npv=500, mu=0, sigma=4, smooth=FALSE, se=FALSE, A=NULL)
{
  if (!is.null(A)){
    tmp=recycle_pv(b, a, first, last, npv=npv, mu=mu, sigma=sigma, A)
    mx = sum(A[last])
    theta = rep(NA,(mx+1))
  }else
  {
    score = possible_scores(a,first,last)
    tmp = pv_recycle(b,a,first,last,score,npv,mu,sigma)
    mx = sum(a[last])
    theta = rep(NA,(mx+1))
  }
  theta[score+1]=rowMeans(tmp)
  if (is.null(A))
  {
    if (smooth) {
      score=0:mx
      theta = predict(lm(theta ~ poly(score,7)))
    }
    sem=rep(NA,(mx+1))
  }else
  {
    if (smooth) {
      score=0:mx
      theta = predict(lm(theta ~ poly(score,7)))
    }
    sem=rep(NA,(mx+1))
  }
  if (se) sem=apply(tmp,1,sd)
  return(list(theta=theta, se=sem))
}

## EAP using Jeffrey's prior: aka propto sqrt(information)
# Uses a weighted average to integrate over a grid defined by:
# grid_from, grid_to and grid_length.
theta_jEAP=function(b, a, first,last, se=FALSE, grid_from=-6, grid_to=6, grid_length=101) 
{
  theta_grid = seq(grid_from, grid_to, length=grid_length)
  f = IJ_(b,a,first,last,theta_grid)
  prior=sqrt(f$I)
  w = sweep(f$l, 2, prior, '*')
  theta = apply(w, 1, function(x)weighted.mean(theta_grid, w=x))
  sem=rep(NA,length(theta))
  if (se)
  {
    f = IJ_(b,a,first,last, theta)
    sem =sqrt((f$I+(f$J/(2*f$I))^2)/f$I^2)
  }
  return(list(theta=theta,se=sem))
}

# computes distribution of score weighted with A conditionally on score weighted with a
# used for estimation of theta from the unweighted score (for instance)
# colSums(G[,,(1-idx)+1]) are elementary symmetric function with a
elsymat <- function(b, a, A, first, last)
{
  ms.a=sum(a[last])
  ms.A=sum(A[last])
  n=length(first)
  G=array(0, c(ms.A+1, ms.a+1, 2))
  # init
  G[1,1,1]=1
  for (j in first[1]:last[1])
  {
    G[A[j]+1,a[j]+1,1]=b[j]
  }
  # recursion
  idx=1
  ms.a=a[last[1]]
  ms.A=A[last[1]]
  for (i in 2:n)
  {
    G[,,idx+1]=G[,,(1-idx)+1]
    for (c in 1:(ms.a+1))
    {
      for (r in 1:(ms.A+1))
      {
        #G[r,c,idx+1]=G[r,c,(1-idx)+1]
        for (j in first[i]:last[i])
        {
          G[r+A[j],c+a[j],idx+1]=G[r+A[j],c+a[j],idx+1]+G[r,c,(1-idx)+1]*b[j]
        }
      }
    }
    idx=1-idx
    ms.a=ms.a+a[last[i]]
    ms.A=ms.A+A[last[i]]
  }
  return(G[,,(1-idx)+1])
}

# ML estimation (using EM) of theta based on A
theta_aA <- function(b,a,A,first,last)
{
  n=length(first)
  ms.a=sum(a[last])
  ms.A=sum(A[last])
  theta=rep(0,ms.A-1)
  G=elsymat(b,a,A,first,last)
  g=colSums(G)
  s_set=which(rowSums(G)>0)-1
  s_set=s_set[-c(1,length(s_set))]
  for (s in s_set)
  {
    converged=FALSE
    while (!converged)
    {
      # E-step (expected score weighted with a, given score weighted with A and current theta)
      p=g*exp((0:ms.a)*theta[s])
      p=p/sum(p)
      E=(G[s+1,]/g)*p
      E=sum((0:ms.a)*E/sum(E))
      # M-step
      converged=(abs(log(E/E_score(theta[s],b,a,first,last)))<1e-6)
      theta[s]=theta[s]+log(E/E_score(theta[s],b,a,first,last))
    }
  }
  #return(c(-Inf,theta,Inf))
  theta = c(-Inf,theta,Inf)
  return(list(theta=theta,se=rep(NA, length(theta))))
}


# Estimate a single ability for a whole score distribution.
# Used for the 3DC standard setting method
# and testing for overdispersion
theta_score_distribution <- function(b,a,first,last,scoretab)
{
  ms.a=sum(a[last])
  theta=0
  np=sum(scoretab)
  escore=-1
  score=(0:ms.a)%*%scoretab
  while (abs(escore-score)>1e-6)
  {
    escore=np*E_score(theta,b,a,first,last)
    theta=theta+log(score/escore)
  }
  return(theta)
}


########################################################
## Score-by-score table. Currently using mean_ElSym
########################################################
# @param m        a rim object produced by fit_inter (but not yet documented anywhere what that looks like)
# @param AB       list: two mutually exclusive subsets of items as indexes first/last/etc.
# @param model    Character: Indicates which model is used: "Rasch" or "IM"
# @return         A list with tbl being a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b
# @details        NA's indicate that a total scores was not possible given the weights
SSTable <- function(m, AB, model) {
  if (model=="IM") {ic=m$est$cIM; b=m$est$bIM} else {ic=m$est$cRM; b=m$est$bRM}
  first = m$inputs$ssI$first
  last =  m$inputs$ssI$last
  C = rep(1:nrow(m$inputs$ssI), m$inputs$ssI$nCat)
  a = m$inputs$ssIS$item_score
  ic = ic[C]
  A = AB[[1]]
  B = AB[[2]]
  ### Check
  if (length(intersect(A,B))!=0) stop("sets not disjunct")

  ## Bookkeeping
  nI=length(first)
  bA=NULL; bB=NULL
  aA=NULL; aB=NULL
  cA=NULL; cB=NULL
  lastA=NULL; firstA=NULL
  lastB=NULL; firstB=NULL
  telAF=1; telBF=1
  for (i in 1:nI)
  {
    if (i %in% A)
    {
      bA=c(bA,b[first[i]:last[i]])
      aA=c(aA,a[first[i]:last[i]])
      cA=c(cA,ic[first[i]:last[i]])
      firstA=c(firstA,telAF)
      lastA=c(lastA,telAF+last[i]-first[i])
      telAF=telAF+last[i]-first[i]+1
    }
    if (i %in% B)
    {
      bB=c(bB,b[first[i]:last[i]])
      aB=c(aB,a[first[i]:last[i]])
      cB=c(cB,ic[first[i]:last[i]])
      firstB=c(firstB,telBF)
      lastB=c(lastB,telBF+last[i]-first[i])
      telBF=telBF+last[i]-first[i]+1
    }
  }
  MscA=sum(aA[lastA])
  MscB=sum(aB[lastB])
  Msc=sum(a[last])
  out=matrix(NA,MscA+1,MscB+1)
  
  ### Rasch Model
  if (model=="RM")
  {
    gA = mean_ElSym(bA,aA,firstA,lastA)
    gB = mean_ElSym(bB,aB,firstB,lastB)
    g =  mean_ElSym(b,a,first,last)
    
    for (s in 0:Msc)
    {
      if (g[s+1]>0)
      {
        for (s_a in max(0,s-MscB):min(s,MscA))
        {
          s_b=s-s_a
          out[s_a+1,s_b+1] = -Inf
          if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
          {
            out[s_a+1,s_b+1] = log(gA[s_a+1]) + log(gB[s_b+1]) - log(g[s+1])
            out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
          }
        }
      }
    }
  }
  
  ### IM
  if (model=="IM")
  {
    logb=log(b); logc=log(ic)
    logbA=log(bA); logcA=log(cA)
    logbB=log(bB); logcB=log(cB)
    for (s in 0:Msc)
    {
      eta=exp(logb+(a*s)*logc)
      etaA=exp(logbA+(aA*s)*logcA)
      etaB=exp(logbB+(aB*s)*logcB)
      gA = mean_ElSym(etaA,aA,firstA,lastA)
      gB = mean_ElSym(etaB,aB,firstB,lastB)
      g =  mean_ElSym(eta,a,first,last)
      if (g[s+1]>0)
      {
        for (s_a in max(0,s-MscB):min(s,MscA))
        {
          s_b=s-s_a
          out[s_a+1,s_b+1] = -Inf
          if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
          {
            out[s_a+1,s_b+1] = log(gA[s_a+1])+log(gB[s_b+1])-log(g[s+1])
            out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
          }
        }
      }
    }
  }
  return(list(tbl=exp(out),m=m,AB=AB,model=model))
}

########################################################
## Score-by-score table ENORM Currently using mean_ElSym
########################################################
# @param m        a parms object produced by fit_enorm
# @param AB       list: two mutually exclusive subsets of items as indexes of items in first and last
# @return         A list with tbl being a score-by-score matrix of probabilities:
#                 P(X^A_+=s_a, X^B_+=s_b|X_+=s) where s=s_a+s_b. 
# @details        NA's indicate that a total scores was not possible given the weights
SSTable_ENORM <- function(b, a, first, last, AB) {
  A = AB[[1]]
  B = AB[[2]]
  ### Check
  if (length(intersect(A,B))!=0) stop("sets not disjunct")
  unionAB = unlist(AB, use.names = F)
  firstA =first[A]; lastA = last[A]
  firstB =first[B]; lastB = last[B]
  MscA = sum(a[lastA])
  MscB = sum(a[lastB])
  Msc  = sum(a[last[unionAB]])
  out = matrix(NA,MscA+1,MscB+1)
  
  gA = mean_ElSym(b,a,firstA,lastA)
  gB = mean_ElSym(b,a,firstB,lastB)
  g =  mean_ElSym(b,a,first[unionAB],last[unionAB])
  
  for (s in 0:Msc)
  {
    if (g[s+1]>0)
    {
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b = s-s_a
        out[s_a+1,s_b+1] = -Inf
        if ((gA[s_a+1]>0)&(gB[s_b+1]>0))
        {
          out[s_a+1,s_b+1] = log(gA[s_a+1]) + log(gB[s_b+1]) - log(g[s+1])
          out[s_a+1,s_b+1] = out[s_a+1,s_b+1] + lchoose(MscA, s_a) + lchoose(MscB, s_b) - lchoose(Msc, s)
        }
      }
    }
  }
  return(list(tbl=exp(out),AB=AB,model="RM"))
}


#################################### Item-Total Regressions


## Wrapper to C function
# using Elsym
ittotmat0 <- function(b,c,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  mm=sum(last-first+1)
  pi=double(mm*(ms+1))
  tmp=.C("ittot_mat0",
         as.double(b),
         as.integer(a),
         as.double(c),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(n),
         as.integer(ms),
         as.double(pi))
  dim(tmp[[8]])=c(mm,(ms+1))
  return(as.matrix(tmp[[8]])) 
}

## low-level check for trivial scores. 
# If TRUE it means that there is no data-reduction
# Each weighted score can be obtained with only one response pattern
# @pi_mat An item-total regression matrix as produced by ittotmat
check_trivial_scores_<-function(pi_mat,scoretab=NULL)
{
  if (!is.null(scoretab)) pi_mat=pi_mat[,scoretab>0]
  pi_mat[(pi_mat==1)|(pi_mat==0)]=NA
  return(all(is.na(pi_mat)))
}

## Wrapper to C function
# currently using mean_ElSym
ittotmat <- function(b,c,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  mm=sum(last-first+1)
  pi=double(mm*(ms+1))
  tmp=.C("ittot_mat",
          as.double(b),
          as.integer(a),
          as.double(c),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(n),
          as.integer(ms),
          as.double(pi))
  dim(tmp[[8]])=c(mm,(ms+1))
  return(as.matrix(tmp[[8]])) 
}



## Get the score distribution of a booklet from fit_enorm
#  If also gives you a smooth version based on a polynomial smoothing
#  of the log-lambda's
#### Example
# db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender=""))
# add_booklet(db, verbAggrData, "agg")
# f=fit_enorm(db)
# xx=dexter:::ENORM2ScoreDist(f,degree=7,booklet_id="data")
# plot(xx$score,xx$p.obs,pch=16, main=paste("degree =", as.character(xx$degree)))
# lines(xx$score,xx$p.obs,col="red")
# lines(xx$score,xx$p.smooth, col="green")
ENORM2ScoreDist <- function (parms, degree=7, booklet_id) 
{
  bk_indx=which(names(parms$inputs$bkList)==booklet_id)
  stopifnot(length(bk_indx)>0) 
  if (parms$inputs$method!="CML") stop("ENORM2ScoreDist only for CML at the moment")
  score_range=0:(length(parms$est$lambda[[bk_indx]])-1)
  a=parms$inputs$ssIS$item_score
  first=parms$inputs$bkList[[bk_indx]]$first
  last=parms$inputs$bkList[[bk_indx]]$last
  lambda=parms$est$lambda[[bk_indx]]
  b=parms$est$b
  log_l=log(lambda)
  
  degree=min(degree,length(!is.na(log_l)))
  qr=lm(log_l~poly(score_range,degree,raw=TRUE)) 
  beta=as.numeric(qr$coefficients)[-1]
  
  lambda[is.na(lambda)]=0
  mx = sum(a[last])
  g = elsym(b,a,first,last)
  sc_obs=vector("numeric", mx+1)
  sc_sm=vector("numeric", mx+1)
  num_obs=0.0
  num_sm=0.0
  for (s in 0:mx)
  {
    sc_obs[s+1]=g[s+1]*lambda[s+1]
    sc_sm[s+1]=g[s+1]*exp(sum(beta*s^(1:degree)))
    num_obs=num_obs+sc_obs[s+1]
    num_sm=num_sm+sc_sm[s+1]
  }
  sc_obs=sc_obs/num_obs
  sc_sm=sc_sm/num_sm
  data.frame(score=score_range,
             n.obs=sc_obs*sum(parms$inputs$stb$N), 
             n.smooth=sc_sm*sum(parms$inputs$stb$N),
             p.obs=sc_obs,
             p.smooth=sc_sm)
}






############################ Elementary Symmetry
# All adapted for inclusion of zero-category
##
elsym <- function(b,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  g=double((ms+1))
  tmp = .C("ElSym",
          as.double(b),
          as.integer(a),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(-1),
          as.integer(-1),
          as.integer(n),
          as.integer(ms),
          as.double(g))
  return(tmp[[9]])
}


#################################### Means esf's
mean_ElSym <- function(b,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  g=double((ms+1))
  tmp = .C("meanElSym",
          as.double(b),
          as.integer(a),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(-1),
          as.integer(-1),
          as.integer(n),
          as.integer(ms),
          as.double(g),
          as.integer(0))
  return(tmp[[9]])
}



