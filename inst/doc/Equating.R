## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(echo = TRUE, message=FALSE)

if (requireNamespace("Cairo", quietly = TRUE)) 
{
   opts_chunk$set(dev='CairoPNG')
}
library(dplyr)
library(ggplot2)

RcppArmadillo::armadillo_throttle_cores(1)

## ----echo=FALSE---------------------------------------------------------------
CurlyBraces <- function(x, y, range, pos = 1, direction = 1 ) {

    a=c(1,2,3,48,50)    # set flexion point for spline
    b=c(0,.2,.28,.7,.8) # set depth for spline flexion point

    curve = spline(a, b, n = 50, method = "natural")$y / 2 

    curve = c(curve,rev(curve))

    a_sequence = rep(x,100)
    b_sequence = seq(y-range/2,y+range/2,length=100)  

    # direction
    if(direction==1)
    a_sequence = a_sequence+curve
    if(direction==2)
    a_sequence = a_sequence-curve

    # pos
    if(pos==1)
    lines(a_sequence,b_sequence) # vertical
    if(pos==2)
    lines(b_sequence,a_sequence) # horizontal

    }

x = seq (-4,4,len=101)
plot (c(1,1-plogis(x,-.5, .5),0), c(1,1-plogis(x,.7,.5),0), type="l", xlim=c(0,1), ylim=c(0,1), xlab='FPR', ylab='TPR',bty='l')
points(1-plogis(0,-.5, .5), 1-plogis(0,.7,.5), pch=19)
points(1-plogis(-1.5,-.5, .5), 1-plogis(-1.5,.7,.5), pch=19, col='gray')

curve(plogis(x,-.5,.5), -4, 4, col=2, xlab='x', ylab='Probability',bty='l')
curve(plogis(x,.7,.5), -4, 4, add=TRUE, col=4)
abline(v=0)
abline(v=-1.5, col='gray')
rg=1-plogis(0,-.5,.5)
text(-1,1-rg/2,'FPR',cex=.6)
CurlyBraces(x=0, y=1-rg/2, range=rg, dir=2)
rg=1-plogis(0,.7,.5)
CurlyBraces(x=0, y=1-rg/2, range=rg, dir=1)
text(1,1-rg/2,'TPR',cex=.6)


## ----echo=TRUE,  results='hide', fig.align='center', fig.height=4, fig.width=4----
library(dexter)

db = start_new_project(verbAggrRules, ":memory:")
add_booklet(db, verbAggrData, "data")

ts = get_testscores(db, item_position < 15) |>
  inner_join(get_testscores(db, item_position >= 15), by='person_id') |>
  rename(ref_test= 'booklet_score.y', new_test = 'booklet_score.x')


ggplot(ts, aes(x = new_test, y = ref_test)) +
  geom_count(show.legend = F) +
  geom_hline(yintercept = 10, colour = 'green') +
  labs(y = 'ref. test score', x = 'target test score') +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

## ----fig.align='center', fig.height=5, fig.width=5----------------------------
prob_pass = tp = rep(0,29)

for (i in seq_along(0:28))
{
  prob_pass[i] = sum(ts$ref_test[ts$new_test==i]>=10) / sum(ts$new_test==i)
  tp[i] = sum(ts$ref_test[ts$new_test>=i]>=10) / sum(ts$new_test>=i)
}

plot(0:28, prob_pass, ylab="Proportion passing the reference test", xlab="New test score", 
     ylim=c(0,1), type = "o", col="red",bty='l')
lines(0:28, tp, type="o", lty=2, col="blue")

## ----fig.align='center', fig.height=5, fig.width=5----------------------------
specificity = sensitivity = rep(0, 29)

for (i in seq_along(0:28))
{
  sensitivity[i] = sum(ts$ref_test[ts$new_test>=i]>=10)/sum(ts$ref_test>=10)
  specificity[i] = sum(ts$ref_test[ts$new_test<i]<10)/(sum(ts$ref_test[ts$new_test<i]<10)+sum(ts$ref_test[ts$new_test>=i]<10)) 
}

plot(0:28, sensitivity, ylab="sensitivity/specificity", xlab="new test score", 
     ylim=c(0,1), type = "o", col="red",bty='l')
lines(0:28, specificity, col="green", type="o")

## ----fig.align='center', fig.height=5, fig.width=5----------------------------
plot(1-specificity, sensitivity, col="green", xlim=c(0,1), ylim=c(0,1), type="l",bty='l')
text(1-specificity, sensitivity, as.character(0:28), cex=0.7, offset = 0)
abline(0,1,lty=2, col="grey")

## ----include=FALSE------------------------------------------------------------
close_project(db)


## -----------------------------------------------------------------------------
#simulate data
n_persons = 700
n_items = 60

theta = c(rnorm(500, 0,2), rnorm(n_persons-500, 0.5,2))
items = data.frame(item_id = sprintf('%02i',1:n_items), item_score=1, beta = runif(n_items,-1,1))

data = r_score(items)(theta)

# incomplete design
data[1:500, 41:60] = NA
data[501:n_persons, 1:20] = NA

## ----fig.align='center', results='hide',fig.height=4,fig.width=4--------------
ref_items = items$item_id[1:40]
target_items = items$item_id[21:60]

p = fit_enorm(data, method='Bayes', nDraws = 5000)

pp = probability_to_pass(data, p,
                         ref_items = ref_items, pass_fail = 23,
                         target_booklets = data.frame(item_id=target_items))
plot(pp, what="equating")

## ----echo=TRUE, fig.align='center', results='hide', fig.height=4, fig.width=4----
plot(pp, what='sens/spec')

## ----include=FALSE------------------------------------------------------------
score_new = coef(pp) |>
  mutate(loss = (1-sensitivity)^2 + (1-specificity)^2) |>
  filter(loss == min(loss)) |>
  pull(score_new)

## ----include=FALSE------------------------------------------------------------
RcppArmadillo::armadillo_reset_cores()

