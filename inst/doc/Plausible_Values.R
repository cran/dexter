## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE)

if (requireNamespace("Cairo", quietly = TRUE)) 
{
  opts_chunk$set(dev='CairoPNG')
}
par_hook = function(before, options, envir)
{
  if(before)
  {
    do.call(par, options$par)
  }
}
knit_hooks$set(par = par_hook)

RcppArmadillo::armadillo_throttle_cores(1)

library(dexter)
library(dplyr)
library(tidyr)

select = dplyr::select # fixes a weird bug

## -----------------------------------------------------------------------------
sim_PV = function(theta, delta) {
  nit = length(delta)
  sim_func = r_score( tibble(item_id=1:nit, item_score=1, delta=delta))
  simulated_data = sim_func(theta)
  parms = fit_enorm(simulated_data)
  pv = plausible_values(simulated_data, parms, nPV = 10)
  
  plot(x=c(min(theta)-.5, max(theta)+.5), y = 0:1, 
       main = paste(nit, "items"), type='n',
       xlab=expression(theta), ylab = expression(Fn(theta)))
  
  select(pv, starts_with('PV')) |>
    lapply(function(x) lines(ecdf(x-mean(x)), col=rgb(.7,.7,.7,.5)))

  lines(ecdf(theta-mean(theta)))
}

## ----fig.align='center', fig.height=4,fig.width=4,par=list(bty='l')-----------
nP = 600
theta = rnorm(nP)
delta = runif(50, -2, 2)
sim_PV(theta, delta)

## ----fig.align='center',  fig.width=7,par=list(bty='l')-----------------------
grp = sample(2, nP, replace = TRUE, prob = c(.5,.5))
theta = rnorm(nP, mean = c(-2,2)[grp], sd = c(0.5,1.5)[grp])
plot(density(theta),main='')

## ----fig.align='center',  fig.width=7,par=list(bty='l',mfrow=c(1,3))----------
sim_PV(theta, delta[1:10])
sim_PV(theta, delta[1:20])
sim_PV(theta, delta)

## -----------------------------------------------------------------------------
sim_PV2 = function(theta, delta) {
  nit = length(delta)
  sim_func = r_score( tibble(item_id=1:nit, item_score=1, delta=delta))
  simulated_data = sim_func(theta)
  parms = fit_enorm(simulated_data, method="Bayes", nDraws = 50)
  
  plot(x=c(min(theta)-.5, max(theta)+.5), y=0:1, main=paste(nit, "items"), type='n',
       xlab=expression(theta), ylab = expression(Fn(theta)))
  
  which.draw = 5*(1:10)
  for (iter in 1:10) {
    pv = plausible_values(simulated_data, parms, parms_draw=which.draw[iter])
    lines(ecdf(pv$PV1-mean(pv$PV1)), col=rgb(.7,.7,.7,.5))
  }
  lines(ecdf(theta-mean(theta)))
}

## ----make_pv2, fig.align='center',fig.width=7, par=list(mfrow=c(1,3),bty='l')----
sim_PV2(theta, delta[1:10])
sim_PV2(theta, delta[1:20])
sim_PV2(theta, delta)

## ----make_pv3, fig.align='center',fig.width=7,par=list(mfrow=c(1,3))----------
sim_PV3 = function(theta, delta) {
  nit = length(delta)
  sim_func = r_score( tibble(item_id=1:nit, item_score=1, delta=delta))
  simulated_data = sim_func(theta)
  parms = fit_enorm(simulated_data)
  pv = plausible_values(simulated_data, parms, nPV = 10, prior_dist = "mixture")
  
  plot(x=c(min(theta)-.5, max(theta)+.5), y=0:1, 
       main=paste(nit, "items"), type='n',
       xlab=expression(theta), ylab = expression(Fn(theta)))
  
  select(pv, starts_with('PV')) |>
    lapply(function(x) lines(ecdf(x-mean(x)), col=rgb(.7,.7,.7,.5)))

  lines(ecdf(theta-mean(theta)))
}

sim_PV3(theta, delta[1:10])
sim_PV3(theta, delta[1:20])
sim_PV3(theta, delta)

## ----make_pv4, fig.align='center',fig.width=7,par=list(mfrow=c(1,3),bty='l')----

sim_PV4 = function(theta, delta, group) {
  nit = length(delta)
  sim_func = r_score( tibble(item_id=1:nit, item_score=1, delta=delta))
  
  # turn simulated data into long format
  simulated_data = sim_func(theta) |>
    get_responses() |>
    inner_join(tibble(person_id=1:length(group), group = group), 
               by='person_id')
  
  parms = fit_enorm(simulated_data)
  
  pv = plausible_values(simulated_data, parms, covariates='group',nPV = 10)
  
  plot(x=c(min(theta)-.5, max(theta)+.5), y=0:1, main=paste(nit, "items"), type='n',
       xlab=expression(theta), ylab = expression(Fn(theta)))
  select(pv, starts_with('PV')) |>
    lapply(function(x) lines(ecdf(x-mean(x)), col=rgb(.7,.7,.7,.5)))
  lines(ecdf(theta-mean(theta)))
}

sim_PV4(theta, delta[1:10],grp)
sim_PV4(theta, delta[1:20],grp)
sim_PV4(theta, delta,grp)


## ----include=FALSE------------------------------------------------------------
RcppArmadillo::armadillo_reset_cores()

