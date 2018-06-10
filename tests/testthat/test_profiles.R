library(dplyr)
library(tibble)
library(DBI)
library(RSQLite)

verbAggCopy = function(pth = '../verbAggression.db')
{
  con = dbConnect(SQLite(), ":memory:")
  db = open_project(pth)
  
  sqliteCopyDatabase(db, con)
  
  dbDisconnect(db)
  return(con)
}


context('Check profile analysis')

test_that('profile analysis verb agg',{
  db = verbAggCopy()
  
  f = fit_enorm(db)
  p = profiles(db, f, 'behavior')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), 0.6, 
            'expected score should have a relation with observed score')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), cor(p$sumScore,p$expected_domain_score),
            'domain should add extra information')
  
  expect_true(all(p %>% 
                    group_by(person_id) %>% 
                    summarise(sum_dif = abs(sum(expected_domain_score) - first(sumScore))) %>%
                    ungroup() %>%
                    pull(sum_dif) < 1e-10), 
              'expected domains scores need to sum to total test score')
  # takes a little long for cran
  skip_on_cran()
  f = fit_enorm(db, method='Bayes')
  p = profiles(db, f, 'behavior')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), 0.6, 
            'expected score should have a relation with observed score (Bayes)')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), cor(p$sumScore,p$expected_domain_score),
            'domain should add extra information (Bayes)')
  
  expect_true(all(p %>% 
                    group_by(person_id) %>% 
                    summarise(sum_dif = abs(sum(expected_domain_score) - first(sumScore))) %>%
                    ungroup() %>%
                    pull(sum_dif) < 1e-10), 
              'expected domains scores need to sum to total test score (Bayes)')
  
  
})



