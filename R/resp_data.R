

# faster factor, does not handle or check for NA values
ffactor = function (x, levels=NULL, as_int=FALSE) 
{
  if(is.null(levels))
  {
    fast_factor(x, as_int)
  } else
  {
    fast_factor_lev(x, levels, as_int)
  }
}



#' Functions for developers
#'
#' These functions are meant for people who want to develop their own models based
#' on the data management structure of dexter. Very little input checking is performed,
#' the benefit is some extra speed over using `get_responses`.
#' Regular users are advised not to use these functions 
#' as incorrect use can easily crash your R-session or lead to unexpected results.
#' 
#'
#' @param dataSrc data.frame, integer matrix, dexter database or `dx_resp_data` object
#' @param qtpredicate quoted predicate
#' @param extra_columns to be returned in addition to person_id, booklet_id, item_score, item_id
#' @param summarised if TRUE, no item scores are returned, just sumscores
#' @param env environment for evaluation of qtpredicate, defaults to caller environment
#' @param protect_x best set TRUE (default)
#' @param retain_person_id whether to retain the original person_id levels or just use arbitrary integers
#' @param merge_within_persons merge different booklets for the same person together
#' @param parms_check data.frame of item_id, item_score to check for coverage of data
#' 
#' @return
#' list of type `dx_resp_data` containing:
#' \describe{
#' \item{when summarised=FALSE}{
#'   x: tibble(person_id, booklet_id, item_id, item_score, booklet_score <, extra_columns>)}
#' \item{when summarised=TRUE}{
#'   x: tibble(person_id, booklet_id, booklet_score <, extra_columns>)}}
#'
#' design: tibble(booklet_id, item_id), sorted
#'
get_resp_data = function(dataSrc, qtpredicate=NULL, 
                         extra_columns=NULL, 
                         summarised=FALSE, env=NULL,
                         protect_x=TRUE, retain_person_id=TRUE,
                         merge_within_persons = FALSE,
                         parms_check=NULL)
{

  if(inherits(dataSrc,'dx_resp_data'))
  {
    if(!is.null(qtpredicate)) 
      stop("predicates don't work in combination with dx_resp_data object yet")
    
    return(resp_data.from_resp_data(dataSrc, extra_columns=extra_columns, summarised=summarised, 
                                    protect_x=protect_x,merge_within_persons = merge_within_persons))
  }
  
  if(is.matrix(dataSrc))
  {
    if(!is.null(qtpredicate)) 
      stop("predicates are not supported for dataSrc of class 'matrix'")
    
    if(!is.null(extra_columns))
      stop("a dataSrc of class 'matrix' is not supported for this function")
    
    return(resp_data.from_matrix(dataSrc, summarised=summarised, retain_person_id=retain_person_id,
                                 merge_within_persons=merge_within_persons, parms_check=parms_check ))
    
  }
  if(is_grouped_df(dataSrc))
    dataSrc = ungroup(dataSrc)
  
  if(is.null(env)) 
    env = caller_env()
  
  if(is_db(dataSrc))
    protect_x = FALSE
  
  booklet_safe = is_bkl_safe(dataSrc, qtpredicate, env) 
  

  # special case that can be done much faster
  # merge_within person could be accomodated, seems not much use though if we lose booklet id
  if(summarised && is_db(dataSrc) && booklet_safe && !merge_within_persons)
  {
    columns = union(extra_columns, c('person_id','booklet_id'))
    
    x = db_get_testscores(dataSrc, qtpredicate=qtpredicate, columns=columns, env=env)
    if(NROW(x)==0)
      stop("no data to analyse")
    x$person_id = ffactor(x$person_id, as_int = !retain_person_id)
    
    if(is.null(qtpredicate))
    {
      design = dbReadTable(dataSrc, 'dxbooklet_design')
      design$booklet_id = ffactor(design$booklet_id)
      x$booklet_id = ffactor(x$booklet_id, levels = levels(design$booklet_id))
    } else
    {
      x$booklet_id = ffactor(x$booklet_id)
      design = dbGetQuery(dataSrc, paste("SELECT booklet_id, item_id, item_position FROM dxBooklet_design",
                                         "WHERE booklet_id IN(",
                                         paste(sql_quote(levels(x$booklet_id),"'"),collapse=','),");"))
      
      design$booklet_id = ffactor(design$booklet_id, levels=levels(x$booklet_id))
    }
    
    design$item_id = ffactor(design$item_id)
    
    if(!is.null(parms_check))
    {
      # is inputs$ssIS a superset of the items and scores that were used to generate the parms?
      
      # parms can be generated in weird unkown ways
      # we know scoring rules is always a superset of responses in the db, so if parms => rules, we're ok
      
      scores = dbGetQuery(dataSrc, "SELECT DISTINCT item_id, item_score FROM dxScoring_rules WHERE item_score>0;") %>%
        semi_join(tibble(item_id = levels(design$item_id)), by='item_id')
      
      # suppress factor warnings, is scores covered by parms (parms is superset)
      # this check covers most cases and should be fast
      suppressWarnings({no_par = dplyr::setdiff(scores, parms_check)})
      
      if(nrow(no_par) > 0)
      {
        # if not all scores are in parms, parms can still be ok if these scores are not in the data either
        # (e.g. because a subset of persons or not all scores in the rules occur in the data)
        # this check is more time consuming but should occur only rarely

        items_in = paste(sql_quote(unique(no_par$item_id), "'"), collapse=',')
        if(is.null(qtpredicate))
        {
          itm_sc = dbGetQuery(dataSrc,
              paste("SELECT DISTINCT item_id, item_score FROM dxResponses INNER JOIN dxScoring_rules USING(item_id,response)",
                          "WHERE item_id IN(",items_in,");"))
        } else
        {
		      pred_sql = qtpredicate_to_sql(qtpredicate, dataSrc, env)	
          itm_sc = try(dbGetQuery(dataSrc,
                          paste("SELECT DISTINCT item_id, item_score FROM", 
                          get_cte(dataSrc, union(c('item_id','item_score'), pred_sql$db_vars)),
                          pred_sql$where,
                          "AND item_id IN(",items_in,");")), silent=TRUE)
          if(inherits(itm_sc,'try-error'))
          {
            itm_sc = dbGetQuery(dataSrc,
                       paste("SELECT DISTINCT ",
                             paste(union(c('item_id','item_score'), pred_sql$db_vars),collapse=','),
                             "FROM", 
                             get_cte(dataSrc, union(c('item_id','item_score'), pred_sql$db_vars)),
                             "WHERE item_id IN(",items_in,");"))
            itm_sc = itm_sc[eval_tidy(qtpredicate, data=itm_sc, env=env) , c('item_id','item_score')]
          }
        }
        
        no_par = dplyr::intersect(itm_sc, no_par)

        
        if(nrow(no_par) > 0)
        {
          message("no parameters for these items and/or scores")
          print(no_par)
          stop("no parameters for some items and/or scores in your data")
        }
      }
    }
    out = list(x=x, design=design, summarised=TRUE, merge_within_persons = merge_within_persons)
    
    class(out) = append('dx_resp_data', class(out))
    return(out)
  }
  
  
  
  columns = unique(c(extra_columns, 'person_id','booklet_id','item_id','item_score'))
  
  if(inherits(dataSrc,'data.frame'))
    columns = intersect(colnames(dataSrc), columns)
  
  x = get_responses_(dataSrc, qtpredicate = qtpredicate, columns = columns, env = env) 
  
  if(NROW(x)==0)
    stop("no data to analyse")
  
  if(inherits(dataSrc,'data.frame'))
  {
    if(nrow(dataSrc) > nrow(x))
      protect_x = FALSE
    mn = min(x$item_score)
    if(is.na(mn))
      stop('item scores may not be NA')
    if(mn<0)
      stop('all item scores must be >= 0')
  }

  
  # very common case, saves .5 seconds with letting from_df sort it out
  if(booklet_safe && is_db(dataSrc) && !summarised)
  {
    x$person_id = ffactor(x$person_id, 
                          dbGetQuery(dataSrc,'SELECT person_id FROM dxpersons 
                                              ORDER BY person_id;')$person_id, 
                          as_int=!retain_person_id)
    
    design = dbGetQuery(dataSrc, "SELECT booklet_id, item_id, item_position FROM dxbooklet_design;")
    
    if(is.null(qtpredicate))
    {
      # assume there are no empty booklets, if there are, sufstats, etc. might be a little less efficient
      design$booklet_id = ffactor(design$booklet_id)
      x$booklet_id = ffactor(x$booklet_id, levels(design$booklet_id))
    } else
    {
      x$booklet_id = ffactor(x$booklet_id)
      design$booklet_id = ffactor(design$booklet_id, levels(x$booklet_id))
      design = filter(design, !is.na(.data$booklet_id))
    }
    
    design$item_id = ffactor(design$item_id)
    x$item_id = ffactor(x$item_id, levels(design$item_id))
    
    
    if(!is.null(parms_check))
    {
      suppressWarnings({
        itm_scr = dbGetQuery(dataSrc,"SELECT DISTINCT item_id, item_score FROM dxscoring_rules WHERE item_score > 0;") %>%
          semi_join(tibble(item_id=levels(design$item_id)), by='item_id') %>%
          anti_join(parms_check, by=c('item_id','item_score'))
      })
      if(NROW(itm_scr)>0)
      {
        itm_scr$item_id = ffactor(itm_scr$item_id, levels=design$item_id)
        itm_scr = itm_scr %>%
          anti_join(x, by=c('item_id','item_score'))
      }
      if(NROW(itm_scr)>0)
      {
        message("no parameters for these items and/or scores")
        print(itm_scr)
        stop("no parameters for some items and/or scores in your data")
      }
    }
    
    if(!is_person_booklet_sorted(x$booklet_id, x$person_id))
    {
      protect_x = FALSE
      x = arrange(x, .data$person_id, .data$booklet_id) 
    }

    if(merge_within_persons)
    {
      # it is likely we have to merge stuff
      # wrapping in function might make R copy more stuff again
      bmap = merge_booklets(x$booklet_id, x$person_id, 
                              design$booklet_id, nlevels(x$booklet_id))
      
      design = design %>%
        rename(old_booklet_id = "booklet_id") %>%
        inner_join(bmap, by='old_booklet_id')
      
      nr = nrow(design)
      
      design = design %>%
        distinct(.data$booklet_id, .data$item_id)
      
      if(nr > nrow(design))
        stop("at least one person has answered at least one item more than once, this is not allowed")
      
      lvls = bmap %>%
        group_by(.data$booklet_id) %>%
        summarize(lev = paste0('(',paste0(.data$old_booklet_id, collapse=', '),')')) %>%
        arrange(.data$booklet_id) %>%
        pull('lev')
      
      attr(x$booklet_id,'levels') = lvls
      class(design$booklet_id) = 'factor'
      levels(design$booklet_id) = lvls
    }
    
    x$booklet_score = mutate_booklet_score(x$person_id, x$booklet_id, x$item_score)
    
    out = list(x=x, design=design, summarised=FALSE, merge_within_persons = merge_within_persons)
    class(out) = append('dx_resp_data', class(out))
    return(out)
  }
  

  return(
    resp_data.from_df(x, extra_columns=extra_columns, 
                      summarised = summarised,
                      booklet_safe = booklet_safe, protect_x = protect_x,
                      merge_within_persons = merge_within_persons, retain_person_id=retain_person_id,
                      parms_check=parms_check))
  
}


resp_data.from_resp_data = function(rsp, extra_columns=NULL, summarised=FALSE, protect_x=TRUE,
                                    merge_within_persons = merge_within_persons) 
{
  
  if(merge_within_persons && !rsp$merge_within_persons)
  {
    return(
      resp_data.from_df(rsp$x, extra_columns=extra_columns, 
                        summarised = summarised,
                        booklet_safe = TRUE, protect_x = protect_x,
                        merge_within_persons = merge_within_persons))
    
  }
  
  # summarising
  if(summarised == rsp$summarised)
    return(rsp)
  
  if(rsp$summarised && !summarised)
    stop("cannot unscramble an egg")
  
  if(protect_x)
  {
    rsp$x$person_id = duplicate(rsp$x$person_id)
    rsp$x$booklet_id = duplicate(rsp$x$booklet_id)
    rsp$x$item_id = duplicate(rsp$x$item_id)
    rsp$x$item_score = duplicate(rsp$x$item_score)
  }
   
  np = summarise_booklet_score(rsp$x$person_id, rsp$x$booklet_id, rsp$x$item_id, rsp$x$item_score)
  if(!is.null(extra_columns))
  {
    #item_id has become permutation index
    indx = head(rsp$x$item_id, np)
    rsp$x[1:np, extra_columns] = rsp$x[indx, extra_columns]
  }
  rsp$x = head(rsp$x[,union(c('booklet_id','person_id','item_score'), extra_columns)],np) %>%
    rename(booklet_score='item_score')
  
  
  rsp$summarised = TRUE
  
  return(rsp)
}


# data.frame x
# set protect_x to FALSE for a possible small speed gain if x may be destroyed
resp_data.from_df = function(x, extra_columns=NULL, summarised=FALSE,
                             booklet_safe = FALSE, protect_x = FALSE,
                             merge_within_persons = FALSE, retain_person_id=TRUE,
                             parms_check=NULL) 
{
  if(NROW(x)==0)
    stop("no data to analyse")
  
  if(!all(c('person_id','item_id','item_score') %in% colnames(x)))
    stop("data should contain the columns 'item_id','item_score','person_id'")
  
  all_columns = intersect(c('person_id','item_id','item_score','booklet_id','item_position', 
                            extra_columns),
                          colnames(x))
  x = x[,all_columns]
  
  pointers = lapply(x, ppoint)
  
  if(!is.factor(x$item_id))
    x$item_id = ffactor(x$item_id)
  
  x$item_score = as.integer(x$item_score) 
  
  if(!is.null(parms_check))
  {
    suppressWarnings({
      uncalibrated = dplyr::setdiff(x[,c('item_id','item_score')], 
                                    parms_check[,c('item_id','item_score')])})
    
    if(nrow(uncalibrated) > 0)
    {
      message("no parameters for these items and/or scores:")
      print(uncalibrated)
      stop("no parameters for some items and/or scores in your data")
    }
  }
  
  if(!is.integer(x$person_id) && !is.factor(x$person_id))
    x$person_id = ffactor(x$person_id, as_int = !retain_person_id)

  if('booklet_id' %in% all_columns)
  {
    if(!is.factor(x$booklet_id))
      x$booklet_id = ffactor(x$booklet_id)
    
    if(merge_within_persons)
    { 
      if(is.unsorted(x$person_id))
        x = arrange(x, .data$person_id) 
    }else if(!is_person_booklet_sorted(x$booklet_id, x$person_id))
    {
      x = arrange(x, .data$person_id, .data$booklet_id)
    }
  } else
  {
    if(is.unsorted(x$person_id))
      x = arrange(x, .data$person_id) 

    x$booklet_id = 1L
    class(x$booklet_id) = 'factor'
    levels(x$booklet_id) = 'bkl'
    booklet_safe = FALSE
  }
  
  if(summarised && protect_x)
  {
    if(ppoint(x$item_score) != pointers$item_score)
      x$item_score = duplicate(x$item_score)
    
    if(ppoint(x$person_id) != pointers$person_id)
      x$person_id = duplicate(x$person_id)
  }  
  
  
  if(booklet_safe && !merge_within_persons)
  {
    design = get_design_C(x$booklet_id, x$item_id)
    if(summarised)
    {
      np = summarise_booklet_score(x$person_id, x$booklet_id, x$item_id, x$item_score)
      if(!is.null(extra_columns))
      {
        #item_id has become permutation index
        indx = head(x$item_id, np)
        x[seq_len(np), extra_columns] = x[indx, extra_columns]
      }
      names(x)[names(x)=='item_score'] = 'booklet_score'
      x = x[seq_len(np),union(c('booklet_id','person_id','booklet_score'), extra_columns)]
      
    } else
    {
      x$booklet_score = mutate_booklet_score(x$person_id, x$booklet_id, x$item_score)
    }
    
  } else
  {
    if(protect_x && ppoint(x$item_score) != pointers$item_score)
      x$item_score = duplicate(x$item_score)
    
    if(summarised)
    {
      res = make_booklets_summed(x$person_id, x$booklet_id, x$item_id, x$item_score, merge_within_persons)
      names(x)[names(x)=='item_score'] = 'booklet_score'
      if(!is.null(extra_columns))
      {
        #item_id has become permutation index
        x[seq_len(res$np), extra_columns] = x[x$item_id[seq_len(res$np)], extra_columns]
      }
      x = x[seq_len(res$np), union(c('booklet_id','person_id','booklet_score'), extra_columns)]
    } else
    {
      x$booklet_score = 0L
      res = make_booklets(x$person_id, x$item_id, x$item_score, x$booklet_id, x$booklet_score, merge_within_persons)
    }
    
    design = res$design
    if('booklet_id' %in% all_columns && !merge_within_persons)
    {
      lvls = res$map_booklet %>%
        arrange( .data$booklet_id) %>%
        pull(.data$org_booklet_id) %>%
        as.character()

      if(anyDuplicated(lvls))
      {
        # one to many, e.g. response!='NA'
        spr = paste0("%0",ceiling(log10(length(lvls))),'i-%s')
        lvls = sprintf(spr,1:length(lvls), lvls)
      } 
    } else
    {
      lvls = sprintf(paste0("bk%0",ceiling(log10(nrow(res$map_booklet))),'i'), 1:nrow(res$map_booklet))
    }
    class(x$booklet_id) = 'factor'
    class(design$booklet_id) = 'factor'
    levels(x$booklet_id) = lvls
    levels(design$booklet_id) = lvls
  }   
  
  
  out = list(x = x,  design = design, summarised = summarised, merge_within_persons = merge_within_persons)
  class(out) = append('dx_resp_data', class(out))
  
  out
}


resp_data.from_matrix = function(X, summarised = summarised, retain_person_id = retain_person_id,
                                 merge_within_persons = merge_within_persons, parms_check = parms_check )
{

  if(merge_within_persons)
    stop('merge within person not yet implemented for matrix dataSrc')
  
  if(!is.numeric(X))
    stop("dataSrc must be a numeric matrix")
  
  if(typeof(X) == 'double')
    mode(X) = 'integer'

  
  if(is.null(colnames(X)))
  {
    colnames(X) = sprintf('item%04i',1:ncol(X))
  } else if(anyDuplicated(colnames(X)))
  {
    stop('dataSrc must not have duplicate column names')
  }
  
  
  if(min(X, na.rm=TRUE) < 0)
    stop("item_scores must be positive numbers")
  
  if(!is.null(parms_check))
  {
    maxs = max(X, na.rm=TRUE) 
    # item_score <= maxs is necessary to prevent out of bounds in C
    parms_check = parms_check %>%
      mutate(item_id = ffactor(as.character(.data$item_id), levels = sort(colnames(X)),as_int=TRUE)) %>%
      filter(!is.na(.data$item_id) & .data$item_id <= ncol(X) & .data$item_score <= maxs) %>%
      arrange(.data$item_id)
    
    if(n_distinct(parms_check$item_id) < ncol(X))
      stop("dataSrc contains items that are not present in the item parameters")
    
    if(!parms_is_superset_matrix(X, parms_check$item_id, parms_check$item_score, maxs))
      stop("your data contains scores that are not present in the item parameters")
  }
  
  
  if(summarised)
  {
    out = make_booklets_summed_matrix(X, ncol(X), nrow(X))
    out$x$person_id = 1:nrow(X)
  } else
  {
    out = make_booklets_matrix(X, ncol(X), nrow(X))

    class(out$x$item_id) = 'factor' 
    levels(out$x$item_id) = colnames(X)
  }
  
  nbk = max(out$design$booklet_id)
  bkstr = paste0('bk%0', ceiling(pmin(log10(nbk),3)),'i')
  
  class(out$x$booklet_id) = class(out$design$booklet_id) = 'factor' 
  levels(out$x$booklet_id) = levels(out$design$booklet_id) = sprintf(bkstr,1:nbk)
  
  
  class(out$design$item_id) =  'factor' 
  levels(out$design$item_id) = colnames(X)
  
  if(retain_person_id && !is.null(rownames(X)))
  {
    class(out$x$person_id) = 'factor'
    levels(out$x$person_id) = rownames(X)
  } 
  out$summarised = summarised
  class(out) = append('dx_resp_data', class(out))
  out
}


# ___________________ methods ______________________ #

intersection = function(respData)
{
  
  if(respData$summarised)
    stop("intersection on summarised data is impossible")
  
  nb = n_distinct(respData$design$booklet_id)
  if(nb>1)
  {
    items = respData$design %>%
      count(.data$item_id) %>%
      filter(.data$n == nb)
    
    if(nrow(items) == 0)
      stop('There are no items in your data that were answered by all respondents')
    
    respData$design = select(items, -.data$n)
    
    respData$design$booklet_id = 1L
    class(respData$design$booklet_id) = "factor"
    attr(respData$design$booklet_id, 'levels') = "intersection"
    
    respData$x = semi_join(respData$x, respData$design, by='item_id')
    
    respData$x$booklet_id = 1L
    class(respData$x$booklet_id) = "factor"
    attr(respData$x$booklet_id, 'levels') = "intersection"
    
    respData$x$booklet_score = mutate_booklet_score(respData$x$person_id, respData$x$booklet_id, respData$x$item_score)
  }
  respData
}

# item_property must be a string indicating a column in x
# no protection against making meaningless combinations if more than one booklet!
# currently removes any other person or item property columns from x
# to~do: extend to keep person properties for future extended use, accept item property in design or as extra
polytomize = function(respData, item_property, protect_x=TRUE)
{
  if(!is.factor(respData$x[[item_property]]))
  {
    respData$x[[item_property]] = ffactor(respData$x[[item_property]]) # knowing the levels would make this faster
  } else if(protect_x)
  {
    respData$x[[item_property]] = duplicate(respData$x[[item_property]])
  }
  if(protect_x)
  {
    respData$x$booklet_id = duplicate(respData$x$booklet_id)
    respData$x$person_id = duplicate(respData$x$person_id)
    respData$x$item_score = duplicate(respData$x$item_score)
    respData$x$booklet_score = duplicate(respData$x$booklet_score)
  }
  
  res = polytomize_C(respData$x$booklet_id, respData$x$person_id, respData$x[[item_property]], 
                       respData$x$item_score, respData$x$booklet_score,
                     nlevels(respData$x[[item_property]]), nlevels(respData$design$booklet_id))
  
  respData$x = respData$x[seq_len(res$n), c('booklet_id','person_id',item_property,'item_score','booklet_score')]
  
  names(respData$x)[names(respData$x) == item_property] = 'item_id'
  
  respData$design = res$design
  respData
}

# filter a dx_resp_data object similar to a dplyr filter statement
# This is less flexible than get_resp_data because
# the predicate here may only contain columns that appear in both the design and the respons data
# on the plus side, if you already have a dx_resp_data object this is much faster than 
# retrieving data from the database
# @parameter respData an object of type dx_resp_data
# @parameter predicate statement (raw/unquoted) to filter data, the predicate may only contain columns that appear in both 
# the design and the data 
# @parameter .recompute_sumscores shall booklet_scores be recomputed
filter.dx_resp_data = function(respData, predicate, env=NULL, .recompute_sumscores = FALSE)
{
  if(is.null(env)) env = caller_env()
  qtpredicate = eval(substitute(quote(predicate)))
  
  if(.recompute_sumscores & respData$summarised) 
    stop('cannot recompute booklet_scores on summarised data')
  
  respData$filter = c(respData$filter, as.character(qtpredicate))
  
  # this works but not with the .data pronoun unfortunately
  # respData$design = respData$design[with(respData$design, eval(partial_eval(qtpredicate,env=env))),]
  # respData$x = respData$x[with(respData$x, eval(partial_eval(qtpredicate,env=env))),]
  
  respData$design = respData$design[eval_tidy(qtpredicate, data = respData$design, env = env),]
  respData$x = respData$x[eval_tidy(qtpredicate, data = respData$x, env = env),]
  
  
  if(.recompute_sumscores)
    respData$x$booklet_score = mutate_booklet_score(respData$x$person_id, respData$x$booklet_id,
                                                     respData$x$item_score) 

  return(respData)
}

# filter join for a dx_resp_data object similar to a dplyr semi_join statement
# @parameter respData an object of type dx_resp_data
# @parameter y tibble for semi_join
# @parameter by may only contain columns that are present in the design and in the data(x)
# @parameter .recompute_sumscores shall booklet_scores be recomputed
semi_join.dx_resp_data = function(respData, y, by, .recompute_sumscores = FALSE)
{
  if(.recompute_sumscores && respData$summarised) 
    stop('cannot recompute booklet_scores on summarised data')
  
  respData$design = semi_join(respData$design, y, by = by)
  respData$x = semi_join(respData$x, y, by = by)

  if(.recompute_sumscores)
    respData$x$booklet_score = mutate_booklet_score(respData$x$person_id, respData$x$booklet_id, 
                                                    respData$x$item_score)
  
  return(respData)
}

# filter join for a dx_resp_data object similar to a dplyr semi_join statement
# @parameter respData an object of type dx_resp_data
# @parameter y tibble for semi_join
# @parameter by may only contain columns that are present in the design and in the data(x)
# @parameter .recompute_sumscores shall booklet_scores be recomputed
anti_join.dx_resp_data = function(respData, y, by, .recompute_sumscores = FALSE)
{
  if(.recompute_sumscores && respData$summarised) 
    stop('cannot recompute booklet_scores on summarised data')
  
  respData$design = anti_join(respData$design, y, by = by)
  respData$x = anti_join(respData$x, y, by = by)
  
  if(.recompute_sumscores)
    respData$x$booklet_score = mutate_booklet_score(respData$x$person_id, respData$x$booklet_id, 
                                                    respData$x$item_score)

  
  return(respData)
}



# INDICES must be one string that indicates column in respData$x
by.dx_resp_data = function (data, INDICES, FUN, ..., simplify = TRUE) 
{
  dsg = data$design
  smr = data$summarised
  args = list(...)
  by(data$x, pull(data$x, INDICES), function(x)
  {
    r = list(summarised = smr, design = dsg, x=x)
    class(r) = append('dx_resp_data', class(r))
    do.call(FUN, append(list(r),args))
  }, 
  simplify = simplify)
}



# compute new levels for item_id in respData
# returns items in same order with new factor levels (underlying integer is different)
# caller should replace levels of item_id in respData by new levels of items
re_factor_item_id = function(respData, items)
{
  rd_lev = levels(respData$design$item_id)
  itm_lev = levels(items)
  if(length(rd_lev) != length(itm_lev) || !all(rd_lev==itm_lev))
  {
    new_lev = c(rd_lev, setdiff(itm_lev,rd_lev))
    items = match(itm_lev, new_lev)[as.integer(items)]
    class(items) = 'factor'
    levels(items) = new_lev
  }
  items
}