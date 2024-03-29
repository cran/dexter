context('test plausible_scores')

library(dplyr)


RcppArmadillo::armadillo_throttle_cores(1)

test_that('plausible scores works',{
  skip_on_cran()
  # simulate some data
  items = tibble(item_id = paste0('i',1:60), delta = runif(60, -2, 2))
  design = data.frame(booklet_id = sort(rep_len(paste0('b',1:5),120)), item_id = paste0('i',1:60), stringsAsFactors = F)
  persons = data.frame(person_id = 1:5000, booklet_id = paste0('b',1:5), stringsAsFactors = F, theta = rnorm(5000))
  responses = persons |>
    inner_join(design, by='booklet_id', relationship = "many-to-many") |>
    inner_join(items, by ='item_id') |>
    mutate(item_score = as.integer(rlogis(n(), theta-delta) > 0)) |>
    select(person_id, booklet_id, item_id, item_score)
  
  
  f = fit_enorm(responses)

  ps_b2 = plausible_scores(responses, parms = f, items = filter(design, booklet_id == 'b2') )
  
  sc_b2 = persons |>
    mutate(booklet_id = 'b2') |>
    inner_join(design, by='booklet_id',relationship = "many-to-many") |>
    inner_join(items, by ='item_id') |>
    mutate(item_score = as.integer(rlogis(n(), theta-delta) > 0)) |>
    group_by(person_id) |>
    summarise(test_score = sum(item_score)) |>
    inner_join(ps_b2, by='person_id')

  

  expect_true(cor(sc_b2$PS1, sc_b2$test_score)>0.65, 
              info=paste('expected correlation test_score and PS>0.65, found:',cor(sc_b2$PS1, sc_b2$test_score)))
  
  
  # keep.observed should work
  ps = plausible_scores(responses) |> 
    rename(PS_keep_true = 'PS1') |>
    inner_join(plausible_scores(responses, keep.observed = FALSE), by= c('person_id','booklet_id')) |>
    inner_join(responses |> 
                 group_by(person_id, booklet_id) |> 
                 summarise(booklet_score = sum(item_score)), 
               by=c('person_id','booklet_id'))
  
  expect_gt(cor(ps$PS_keep_true, ps$booklet_score), cor(ps$PS1, ps$booklet_score)+.05)
  
  expect_false(any(ps$PS_keep_true<ps$booklet_score))
  
  expect_lt(abs(mean(ps$PS1)-mean(ps$PS_keep_true)),0.15)
  
})

RcppArmadillo::armadillo_reset_cores()