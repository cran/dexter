% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plausible_values.R
\name{plausible_values}
\alias{plausible_values}
\title{Draw plausible values}
\usage{
plausible_values(
  dataSrc,
  parms = NULL,
  predicate = NULL,
  covariates = NULL,
  nPV = 1,
  parms_draw = c("sample", "average"),
  prior_dist = c("normal", "mixture"),
  merge_within_persons = FALSE
)
}
\arguments{
\item{dataSrc}{a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score}

\item{parms}{An object returned by function \code{fit_enorm} containing parameter estimates or a data.frame with columns item_id, item_score and, 
beta. If parms are provided, item parameters are considered known. If parms is NULL, they will be estimated Bayesianly.}

\item{predicate}{an expression to filter data. If missing, the function will use 
all data in dataSrc}

\item{covariates}{name or a vector of names of the variables to group the populations used to improve the prior.
A covariate must be a discrete person property (e.g. not a float) that indicates nominal categories, e.g. gender or school.
If dataSrc is a data.frame, it must contain the covariate.}

\item{nPV}{Number of plausible values to draw per person.}

\item{parms_draw}{when the item parameters are estimated with method "Bayes" (see: \code{\link{fit_enorm}}), 
parms_draw specifies whether to use a sample (a different item parameter draw for each plausible values draw) or the posterior mean
of the item draws. Alternatively, it can be an integer specifying a specific draw. It is ignored when parms is not estimated Bayesianly.}

\item{prior_dist}{use a normal prior for the plausible values or a mixture of two normals. 
A mixture is only possible when there are no covariates.}

\item{merge_within_persons}{If a person took multiple booklets, this indicates
whether plausible values are generated per person (TRUE) or per booklet (FALSE)}
}
\value{
A data.frame with columns booklet_id, person_id, booklet_score, any covariate columns, and nPV plausible values
named PV1...PVn.
}
\description{
Draws plausible values based on test scores
}
\details{
When the item parameters are estimated using \code{fit_enorm(..., method='Bayes')} and parms_draw = 'sample', the uncertainty 
of the item parameters estimates is taken into account when drawing multiple plausible values. 

In there are covariates, the prior distribution is a hierarchical normal with equal variances across groups. When there is only
one group this becomes a regular normal distribution. When there are no covariates and prior_dist = "mixture", the prior is a mixture
distribution of two normal distributions which gives a little more flexibility than a normal prior.
}
\examples{

\dontshow{ RcppArmadillo::armadillo_throttle_cores(1)}

db = start_new_project(verbAggrRules, ":memory:", 
   person_properties=list(gender="<unknown>"))
add_booklet(db, verbAggrData, "agg")
add_item_properties(db, verbAggrProperties)

f=fit_enorm(db)
pv_M=plausible_values(db,f,(mode=="Do")&(gender=="Male"))
pv_F=plausible_values(db,f,(mode=="Do")&(gender=="Female"))

par(mfrow=c(1,2))

plot(ecdf(pv_M$PV1), 
   main="Do: males versus females", xlab="Ability", col="red")
lines(ecdf(pv_F$PV1), col="green")
legend(-2.2,0.9, c("female", "male") , 
   lty=1, col=c('green', 'red'), bty='n', cex=.75)

pv_M=plausible_values(db,f,(mode=="Want")&(gender=="Male"))
pv_F=plausible_values(db,f,(mode=="Want")&(gender=="Female"))

plot(ecdf(pv_M$PV1), 
   main="Want: males versus females", xlab=" Ability", col="red")
lines(ecdf(pv_F$PV1),col="green")
legend(-2.2,0.9, c("female", "male") , 
   lty=1, col=c('green', 'red'), bty='n', cex=.75)
   
   
close_project(db)    

\dontshow{ RcppArmadillo::armadillo_reset_cores()}

}
\references{
Marsman, M., Maris, G., Bechger, T. M., and Glas, C.A.C. (2016) What can we learn from plausible values? 
Psychometrika. 2016; 81: 274-289. See also the vignette.
}
