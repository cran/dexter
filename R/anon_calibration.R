######################################################################
# Estimation of Rasch and Interaction Model 
#######################################################################
# Using Newton-Raphson per item
# Currently with elsym and scale
# @param ss list containing:
#   il: one row per item, ordered item_id
#       tibble(first, last, sufC <sum(item_score*sumScore)>, nCat <nbr of score categories including 0>)
#   sl: one row per item-scorecat(including 0), ordered item_id, item_score
#       tibble(item_score, sufI <count for item_score>, sufC <sum(item_score * sumScore)>)
#   tl: one row per testscore (only one test allowed), complete range 0:max_observed, ordered by test_score
#       tibble(N <count for test score>)
# @returns:
#      bRM:    Parameter estimates of the Rasch model
#      bIM:    Parameter estimates of the Interaction model
#      cIM:    Estimate of (log-)interaction parameter
#      cRM:    Interaction parameters under the Rasch model: all equal to 1
#      HRM:    Block diag. Asymptotic var-covvar matrix of item parameters under RM
#     se.c:   Standard error of interaction parameter
# fit.stat: log(cIM)/se.c. Wald statistic normally distributed under Rasch model
#########################################################################

EstIM  <- function(first,last, nCat, a, sufI, sufC, scoretab) {
  #first = ss$il$first
  #last = ss$il$last
  #a = ss$sl$item_score
  #sufI = ss$sl$sufI
  #sufC = ss$il$sufC
  #C = rep(1:nrow(ss$il), ss$il$nCat)
  #scoretab = ss$tl$N
  C = rep(1:length(first), nCat)
  
  m=sum(scoretab) ##
  nI=length(last)
  b=rep(0,length(sufI))
  ic=rep(1,nI)
  var.ic=vector("numeric", nI)
  HRM=matrix(0,length(b),length(b))
  
  # Identification
  b[sufI>0]=1
  upd_set=vector("list",nI)
  for (i in 1:nI)
  {
    upd_set[[i]]=which(sufI[first[i]:last[i]]>0)
    upd_set[[i]]=upd_set[[i]][-1]
    upd_set[[i]]=(first[i]:last[i])[upd_set[[i]]]
  }
  
  if (check_trivial_scores_(ittotmat0(b,ic[C],a,first,last),scoretab)){
    warning("Only trivial weighted scores observed") }
  
  converged=2
  scale=2
  while(converged>0.01)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last) ##
    pi_mat[is.na(pi_mat)]=0
    for (i in 1:nI)
    {
      if (length(upd_set[[i]])>0)
      {
        pi=pi_mat[upd_set[[i]],,drop=FALSE]
        E=sufI[upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%tcrossprod(diag(scoretab), pi) #diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # NR update for parameters of item i
        update=solve(H*scale,E)
        b[upd_set[[i]]]=b[upd_set[[i]]]*exp(update)
        converged=max(converged,max(abs(E))/m) #
        HRM[upd_set[[i]],upd_set[[i]]]=H
      }
    }
    if (converged<1) scale=1
  }
  
  bRM=b
  cRM=ic
  
  ## IM
  converged=2
  scale=2
  while(converged>0.001)
  {
    converged=-1
    pi_mat=ittotmat0(b,ic[C],a,first,last)
    pi_mat[is.na(pi_mat)]=0
    for (i in 1:nI)
    {
      # gradient and hessian for thresholds of item i
      if (length(upd_set[[i]])>0)
      {
        pi = pi_mat[upd_set[[i]], , drop=FALSE]
        E = sufI[upd_set[[i]]] - pi %*% scoretab
        H = -pi %*% tcrossprod(diag(scoretab), pi)
        diag(H) = pi %*% scoretab + diag(H)
        
        # gradient and hessian for interaction parameter
        ncol_pi = ncol(pi)
        nrow_pi = nrow(pi)
        
        E = c(E, sufC[i])
        H = cbind(H, rep.int(0, nrow(H)))
        H = rbind(H, rep.int(0, ncol(H)))
        k = 1
        e0 = 0; e1 = 0
        f = matrix(0, nrow_pi, ncol_pi)
        g = matrix(0, nrow_pi, ncol_pi)
        h = 0
        for (j in upd_set[[i]])
        {
          E[length(E)]=E[length(E)]-a[j]*sum((0:(ncol_pi-1))*scoretab*pi[k,])
          e0=e0+a[j]*pi[k,]
          e1=e1+a[j]^2*pi[k,]
          f[k,]=a[j]*(0:(ncol_pi-1))*pi[k,]
          g[k,]=pi[k,]
          h=h+a[j]*(0:(ncol_pi-1))*pi[k,]
          k=k+1
        }
        H[nrow(H),nrow(H)]=sum((0:(ncol_pi-1))^2*(e1-e0^2)*scoretab)
        for (k in 1:nrow(f))
        {
          H[k,nrow(H)]=sum((f[k,]-g[k,]*h)*scoretab)
          H[nrow(H),k]=H[k,nrow(H)]
        }
        # NR update for parameters of item i
        update=solve(H*scale,E)
        b[upd_set[[i]]]=b[upd_set[[i]]]*exp(update[-length(update)])
        ic[i]=ic[i]*exp(update[length(update)])
        var.ic[i]=solve(H)[nrow(H),nrow(H)]
        converged=max(converged,max(abs(E))/m)
      }
    }
    if (converged<1) scale=1
  }
  
  return(list(bRM=bRM,cRM=cRM,bIM=b,cIM=ic,se.c=sqrt(var.ic),HRM=HRM, fit.stats=log(ic)/sqrt(var.ic)))
}



##################################################### calibrate incomplete designs: 
# CML
# Bayes
# Note: if you wish to use mean_Elsym. Replace g=elsym(b, a, booklet[[bl]]$first[-ii], booklet[[bl]]$last[-ii]) by
#   g = mean_ElSym(b, a, booklet[[bl]]$first[-ii], booklet[[bl]]$last[-ii])
#   lg1 = length(g) - 1
#   g = g*choose(lg1, 0:lg1)

#####################################################

#### Bayes
calibrate_Bayes = function(itemList, booklet, sufI, b, a, first, last, nIter, fixed_b=NULL, lambda_out=FALSE) 
{
  nb = length(booklet)
  n = length(itemList)
  y = rep(0, length(sufI))
  z = NULL
  bx = matrix(0, nIter, length(b))
  lx=list()
  if (lambda_out)
  {
    length(lx)=nb
    names(lx)=names(booklet)
    for (bl in 1:nb) lx[[bl]]=matrix(0,nIter,sum(a[booklet[[bl]]$last])+1)
  }
  
  pb = txtProgressBar(min=0, max=nIter)
  for (iter in 1:nIter)
  {
    for (bl in 1:nb)
    {
      # data augmentation
      g=elsym(b, a, booklet[[bl]]$first, booklet[[bl]]$last)
      z[bl] = rgamma(1, shape=booklet[[bl]]$m, 
                     rate=sum(g*booklet[[bl]]$lambda))
      # update lambda
      idx = which(g != 0.0)
      booklet[[bl]]$lambda[idx] = rgamma(length(idx), shape=booklet[[bl]]$scoretab[idx]+0.1, 
                                         rate=(g*z[bl])[idx]) # 1.1
      booklet[[bl]]$lambda[-idx] = 0.0
      # scale lambda such that g*lambda~scoretab
      booklet[[bl]]$lambda[idx] = booklet[[bl]]$m*booklet[[bl]]$lambda[idx]/sum(g*booklet[[bl]]$lambda)
      z[bl] = rgamma(1, shape=booklet[[bl]]$m, 
                     rate=sum(g*booklet[[bl]]$lambda))
    }
    for (i in 1:n)
    {
      y[first[i]:last[i]] = 0.0
      for (bl in itemList[[i]])
      {
        ii = which(booklet[[bl]]$first==first[i])
        g = elsym(b, a, booklet[[bl]]$first[-ii], booklet[[bl]]$last[-ii])
        for (j in first[i]:last[i])
        {
          if (a[j] == 0) y[j]=y[j]+z[bl]*sum(g*head(booklet[[bl]]$lambda,-a[last[i]]))
          if ((a[j] != a[last[i]])&(a[j]!=0)) y[j]=y[j]+z[bl]*sum(g*head(tail(booklet[[bl]]$lambda,-a[j]),-(a[last[i]]-a[j])))
          if (a[j] == a[last[i]]) y[j]=y[j]+z[bl]*sum(g*tail(booklet[[bl]]$lambda,-a[j]))
        }
      }
      b[first[i]:last[i]] = rgamma(1+last[i]-first[i],shape=sufI[first[i]:last[i]]+1.1,
                                   rate=y[first[i]:last[i]]) #1.1
    }
    # identify within items
    for (i in 1:n)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    if (is.null(fixed_b))
    {
      f=b[2] 
      b[-first] = b[-first]/(f^(a[-first]/a[2]))
      # Lambda
      for (bl in 1:nb) 
      {
        booklet[[bl]]$lambda = booklet[[bl]]$lambda*f^(0:sum(a[booklet[[bl]]$last]))
        booklet[[bl]]$lambda = booklet[[bl]]$lambda/booklet[[bl]]$lambda[1]
        if (lambda_out) lx[[bl]][iter,]=booklet[[bl]]$lambda
      }
    }else
    {
      fixed_set=which(!is.na(fixed_b))
      b[fixed_set]=fixed_b[fixed_set]
      if (lambda_out) { for (bl in 1:nb) lx[[bl]][iter,]=booklet[[bl]]$lambda }
    }
    b[is.nan(b)] = 1 # deal with items that are not in data
    bx[iter,] = b
    setTxtProgressBar(pb, value=iter)
  }
  close(pb)
  
  report=toOPLM(a,bx, first, last, H=NULL,fixed_b=fixed_b)
  return(list(a=a, b=report$b_renorm,lambda=lx, beta=report$beta))
}

### Estimate lambda from CML
## arguments first, last and scoretab are for specific booklet
est_lambda <- function(b, a, first, last, scoretab)
{
  ifelse(scoretab>0, scoretab/(elsym(b,a,first,last)*sum(scoretab)), NA) 
}  


#####################################################################
### Do CML
## Must be changed to handle the case where 0 categories does not occur
## Or category-scores (a) are not increasing
## If fixed_b is not NULL unfixed c.q. free parameters are NA
## Note that fixed_b contains values in the dexter parametrisation (including 0 category)

calibrate_CML <- function(booklet, sufI, a, first, last, nIter, fixed_b=NULL) {
  nb = length(booklet)
  ni = length(first)
  EsufI = sufI
  max_nr_iter = 30
  
  if (is.null(fixed_b)) # if no fixed parameters
  {
    nn= sum(sufI)
    b = rep(1,length(a))
    ## Implicit Equations  ###
    converged=FALSE
    iter=0
    pb = txtProgressBar(min=0, max=nIter)
    while ((!converged)&(iter<=nIter))
    {
      iter=iter+1
      EsufI=EsufI-EsufI
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
      }
      b = b*sufI/EsufI
      converged=(max(abs(sufI-EsufI))/nn<1e-04)
      setTxtProgressBar(pb, value=iter)
    }
    ie_iter=iter
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
    ### identification ###
    # within items
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    # between items
    ref_cat=2
    b[-first] = b[-first]/(b[ref_cat]^(a[-first]/a[ref_cat]))
    
    
    ###  NR  ###
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1
    while ((!converged)&&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI=EsufI-EsufI
      H=H-H
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
        H     = H     + H.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
      }
      # identify
      for (i in 1:ni)
      {
        H[first[i],first[i]]=1
        EsufI[first[i]]=sufI[first[i]]
      }
      H[ref_cat,]=0
      H[,ref_cat]=0
      H[ref_cat,ref_cat]=1
      EsufI[ref_cat]=sufI[ref_cat]
      b = b*exp(solve(H*scale,sufI-EsufI))
      converged=(max(abs(EsufI-sufI))/nn<1e-10)
      setTxtProgressBar(pb, value=iter)
      if (nr_iter==2) scale=1
    }
    close(pb)
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  }else  ### if fixed parameters
  {
    fixed_set=which(!is.na(fixed_b))
    update_set=which(is.na(fixed_b))
    b=fixed_b
    ni_free=sum(is.na(fixed_b[last]))
    b[update_set]=1
    nn=0
    for (bl in 1:nb) nn=nn+booklet[[bl]]$m
    nn=nn*ni_free
    
    converged=FALSE
    iter=0
    pb = txtProgressBar(min=0, max=nIter)
    while ((!converged)&(iter<=nIter))
    {
      iter=iter+1
      EsufI=EsufI-EsufI
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
      }
      b[update_set] = b[update_set]*sufI[update_set]/EsufI[update_set]
      converged=(max(abs(sufI[update_set]-EsufI[update_set]))/nn<1e-04)
      setTxtProgressBar(pb, value=iter)
    }
    ie_iter=iter
    if (!converged) warning(paste('Implicit Equations not Converged in',as.character(nIter),"iterations"))
    
    for (i in 1:ni)
    {
      range=first[i]:last[i]
      b[range]=b[range]/b[first[i]]
    }
    
    H=matrix(0,length(a),length(a))
    converged=FALSE
    nr_iter=0
    scale=1
    while ((!converged)&(nr_iter<max_nr_iter))
    {
      iter=iter+1
      nr_iter=nr_iter+1
      EsufI=EsufI-EsufI
      H=H-H
      for (bl in 1:nb)
      {
        EsufI = EsufI + E.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab) 
        H     = H     + H.STEP(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
      }
      # identify
      for (i in 1:length(first))
      {
        H[first[i],first[i]]=1
        EsufI[first[i]]=sufI[first[i]]
      }
      H[fixed_set,]=0
      H[,fixed_set]=0
      diag(H)[fixed_set]=1
      EsufI[fixed_set]=sufI[fixed_set]
      b = b*exp(solve(H*scale,sufI-EsufI))
      converged=(max(abs(EsufI[update_set]-sufI[update_set]))/nn<1e-10)
      setTxtProgressBar(pb, value=iter)
      scale=1
    }
    close(pb)
    if (!converged) warning(paste('Newton-Raphson not Converged in',as.character(nr_iter),"iterations"))
  }
  
  lx=list()
  length(lx)=nb
  names(lx)=names(booklet)
  for (bl in 1:nb) lx[[bl]]=est_lambda(b,a,booklet[[bl]]$first,booklet[[bl]]$last,booklet[[bl]]$scoretab)
  
  report = toOPLM(a, b, first, last, H=H, fixed_b=fixed_b)
  return(list(b=report$b_renorm, H=H, beta=report$beta, acov.beta=report$cov.beta, 
              lambda=lx, n_iter=iter, nr_iter = nr_iter, ie_iter=ie_iter))
}

#################################### Wrappers for C- FUnctions for CML
## THis version uses Elsym and is adapted for inclusion of b_0. However
## THere must still be a small error somewhere.
H.STEP <- function(b,a,first,last,nscore)
{
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)^2)
  tmp=.C("H",
         as.double(b),
         as.integer(a),
         as.integer(length(b)),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[9]]
  tmp=as.matrix(tmp)
  dim(tmp)=c(length(b),length(b))
  tmp=tmp+t(tmp)
  diag(tmp)=diag(tmp)/2
  return(tmp)
}

## This version uses Elsym0... 
H.STEP0 = function(b,a,first,last,nscore)
{
  first=first+1
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)^2)
  tmp=.C("H0",
         as.double(b),
         as.integer(a),
         as.integer(length(b)),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[9]]
  tmp=as.matrix(tmp)
  dim(tmp)=c(length(b),length(b))
  tmp=tmp+t(tmp)
  diag(tmp)=diag(tmp)/2
  return(tmp)
}

########################################
# version with Elsym
E.STEP <- function(b,a,first,last,nscore)
{
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)) 
  tmp=.C("E",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[8]]
  return(tmp)
}

# version with Elsym0
E.STEP0 <- function(b,a,first,last,nscore)
{
  first=first
  n=length(first)
  ms=length(nscore)-1
  output = double(length(b)) 
  tmp=.C("E0",
         as.double(b),
         as.integer(a),
         as.integer(first-1),
         as.integer(last-1),
         as.integer(nscore),
         as.integer(n),
         as.integer(ms),
         as.double(output))[[8]]
  return(tmp)
}
