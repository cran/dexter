
#' Draw plausible test scores
#'
#' Draw plausible, i.e. posterior predictive sumscores on a set of items. 
#' 
#' A typical use of this function is to generate plausible scores on
#' a complete item bank when data is collected using an incomplete design
#'
#' @param dataSrc a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score
#' @param predicate an expression to filter data. If missing, the function will use 
#' all data in dataSrc
#' @param parms An object returned by function \code{fit_enorm} and containing
#' parameter estimates. If parms is given the function provides plausible scores conditional on the 
#' item parameters. These are considered known. If \code{parms} is \code{NULL}, Bayesian parameters are calculated from the datasrc
#' @param parms_draw when the item parameters are estimated Bayesianly (see: \code{\link{fit_enorm}}), 
#' parms_draw specifies whether to use a sample(a different item parameter draw for each plausible values draw) or the posterior mean
#' of the item draws. Alternatively, it can be an integer specifying a specific draw. Ignored when parms is not estimated Bayesianly.
#' @param items vector of item_id's, this specifies the itemset to generate the testscores for. If \code{items} is \code{NULL} 
#' all items occurring in \code{dataSrc} are used.
#' @param covariates name or a vector of names of the variables to group the population, used to update the prior.
#' A covariate must be a discrete person covariate (e.g. not a float) that indicates nominal categories, e.g. gender or school
#' If dataSrc is a data.frame, it must contain the covariate.
#' @param keep.observed If responses to one or more of the items have been observed,
#' the user can choose to keep these observations or generate new ones. 
#' @param nPS Number of plausible testscores to generate per person.
#' @param prior_dist use a normal prior for the plausible values or a mixture of two normals. 
#' A mixture is only possible when there are no covariates.
#' @param merge_within_persons If a person took multiple booklets, this indicates
#' whether plausible scores are generated per person (TRUE) or per booklet (FALSE)
#' @return A data.frame with columns booklet_id, person_id, booklet_score and nPS plausible scores
#' named PS1...PSn.
#'  
plausible_scores = function(dataSrc, parms=NULL, predicate=NULL, items=NULL, parms_draw = c('sample','average'),
                            covariates=NULL, nPS=1, prior_dist = c("normal", "mixture"),
                            keep.observed=TRUE,merge_within_persons=FALSE)  
{
  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  check_dataSrc(dataSrc)
  
  pb = get_prog_bar(nsteps=if(is.null(parms)) 130 else 100, 
                    retrieve_data = is_db(dataSrc))
  on.exit({pb$close()})
  
  
  respData = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, extra_columns=covariates, env=env,
                           merge_within_persons=merge_within_persons)
  
  if(is.null(items))
  {
    items = levels(respData$design$item_id)
  } else if(inherits(items,'data.frame'))
  {
    items = as.character(unique(items$item_id))
  } else
  {
    items = as.character(unique(items))
  }
  
  # if there are no params, all of items must be in data
  # if there are params, all of items must be in params
  
  if(is.null(parms) && !all(items %in% levels(respData$design$item_id)))
  {
    stop_("`items` contains item_id's not found in the data, you must either provide parameters reparately or ",
          "specify only items present in your data")
  } else if(!is.null(parms))
  {
    if(inherits(parms,'data.frame')) parms_items = as.character(unique(parms$item_id))
    else parms_items = unique(coef(parms)$item_id)
    
    if(!all(items %in% parms_items))
      stop_("`items` contains item_id's not found in the parameters")
  }
  
  # generate plausible values and params
  res = plausible_values_(respData, parms=parms, covariates=covariates, 
                          nPV=nPS, parms_draw = parms_draw, 
                          prior_dist = prior_dist)
  
  parms = res$parms
  pv = res$pv
  
  items = factor(items,levels=levels(parms$items$item_id))
  
  fl = parms$items |>
    filter(.data$item_id %in% items) |>
    mutate(first = .data$first-1L, last = .data$last-1L)
  
  a = parms$a
  if(is.matrix(parms$b))
  {
    b = t(parms$b)
    bstep = as.integer((ncol(b)-1)/max(nPS-1,1))
  } else
  {
    b = matrix(parms$b,ncol=1)
    bstep = 0L
  }
  
  if(keep.observed && any(respData$design$item_id %in% items))
  {
    # keep track of sumscore on selected items
    
    respData = semi_join_rd(respData, tibble(item_id=items), by='item_id', .recompute_sumscores = TRUE)
    respData = get_resp_data(respData, summarised = TRUE, protect_x = !is_db(dataSrc))
    
    pv = pv |> 
      select(-'booklet_score') |>
      left_join(respData$x,  by=c("person_id", "booklet_id")) |>
      mutate(booklet_score = coalesce(.data$booklet_score, 0L))
    
    pv = lapply(split(pv, pv$booklet_id), function(pvbk)
    {
      bk = pvbk$booklet_id[1]
      
      fl_bk = fl |>
        anti_join(filter(respData$design, .data$booklet_id == bk), by='item_id')
      
      #nothing to augment case
      if(nrow(fl_bk) == 0)
      {
        for(pn in sprintf('PV%i',1:nPS)) pvbk[[pn]] = pvbk$booklet_score
      } else
      {
        b_index = 1L
        for(pn in sprintf('PV%i',1:nPS))
        {
          pvbk[[pn]] = sampleNRM2_test(pvbk[[pn]], b[,b_index], a, fl_bk$first, fl_bk$last)[,1,drop=TRUE] + pvbk$booklet_score
          b_index = b_index + bstep
        }    
      }
      pvbk
    }) |>
      bind_rows()
    
  } else
  {
    b_index = 1L
    
    for(pn in sprintf('PV%i',1:nPS))
    {
      pv[[pn]] = sampleNRM2_test(pv[[pn]], b[,b_index], a, fl$first, fl$last)[,1,drop=TRUE]
      b_index = b_index + bstep
    }
  }
  
  pv |>
    select(-'booklet_score') |>
    rename_with(gsub, pattern='^PV(?=\\d+$)',replacement='PS', perl=TRUE)  |>
    mutate_if(is.factor, as.character) |>
    df_format()
}

