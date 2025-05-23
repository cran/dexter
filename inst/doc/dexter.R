## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE,fig.align='center', fig.width=6, fig.height=5, message = FALSE)

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

## -----------------------------------------------------------------------------
library(dplyr)
library(dexter)
head(verbAggrRules, 10)

## ----include=FALSE------------------------------------------------------------
db = start_new_project(verbAggrRules, ":memory:", person_properties=list(gender="unknown"))

## -----------------------------------------------------------------------------
add_booklet(db, verbAggrData, "agg")

## -----------------------------------------------------------------------------
add_item_properties(db, verbAggrProperties)

## -----------------------------------------------------------------------------
head(verbAggrProperties)

## -----------------------------------------------------------------------------
get_booklets(db)
head(get_items(db))
get_persons(db) |> 
  glimpse()

## -----------------------------------------------------------------------------
tt = tia_tables(db)

## ----echo=FALSE---------------------------------------------------------------
kable(tt$booklets, digits=3)

## ----echo=FALSE---------------------------------------------------------------
kable(tt$items, digits=3)

## -----------------------------------------------------------------------------
distractor_plot(db, 'S1DoShout')

## -----------------------------------------------------------------------------
m = fit_inter(db, booklet_id=='agg')
plot(m, "S1DoScold", show.observed=TRUE)

## -----------------------------------------------------------------------------
plot(m, 'S1DoCurse', summate=FALSE)

## ----par=list(mfrow=c(2,2))---------------------------------------------------
mSit = fit_domains(db, item_property= "situation")
plot(mSit)

## ----results='hide'-----------------------------------------------------------
parms = fit_enorm(db)

## ----results='hide'-----------------------------------------------------------
parms_gibbs = fit_enorm(db, method='Bayes')

## ----echo=FALSE---------------------------------------------------------------
kable(head(coef(parms_gibbs)), digits=3)

## -----------------------------------------------------------------------------
pv = plausible_values(db, parms)
plot(density(pv$PV1), bty='l', main='verbal aggression', xlab='plausible value')

## ----par=list(bty='n', fg='white')--------------------------------------------
pv = merge(pv, get_persons(db))

boxplot(PV1~gender, data=pv, border='black')

## ----fig.width=5, fig.height=5------------------------------------------------
profile_plot(db, item_property='mode', covariate='gender')

## -----------------------------------------------------------------------------
close_project(db)

RcppArmadillo::armadillo_reset_cores()

