#####################
# Some functions to transform user-provided (i.e., fixed) parameter values from one 
# parameterization to the other.
#####################
beta2eta <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$eta = beta2eta_(df.new$item_score, df.new$beta, first, last)
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(parms.df$item_score, parms.df$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          beta = rep(0,length(tmp$a)), eta = rep(0,length(tmp$a)),
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$beta[(tmp$first[i]+1):tmp$last[i]] = parms.df$beta[first[i]:last[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

beta2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$b = eta2b_(beta2eta_(df.new$item_score, df.new$beta, first, last))
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(df.new$item_score, df.new$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          beta = rep(0,length(tmp$a)), b = tmp$b,
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$beta[(tmp$first[i]+1):tmp$last[i]] = parms.df$beta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

eta2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=FALSE)
{
  df.new = parms.df
  df.new$b = eta2b_(df.new$eta)
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(df.new$item_score, df.new$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          eta = rep(0,length(tmp$a)), b = tmp$b,
                          stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

b2b <-function(first, last, parms.df, out.zero=TRUE, in.zero=TRUE)
{
  df.new=parms.df
  if (out.zero!=in.zero)
  {
    if (in.zero) df.new = parms.df[-first,] # zero in but not out
    if (out.zero) # zero out but not in
    {
      tmp = add_zero(parms.df$item_score, parms.df$b, first, last)
      df.new = data.frame(item_id = rep("i",length(tmp$a)), item_score = tmp$a, 
                          b = tmp$b, stringsAsFactors = FALSE)
      for (i in 1:length(tmp$first))
      {
        df.new$item_id[tmp$first[i]:tmp$last[i]] = parms.df$item_id[first[i]]
        df.new$eta[(tmp$first[i]+1):tmp$last[i]] = parms.df$eta[first[i]:last[i]]
      }
    }
  }
  return(df.new)
}

# item_id, item_score, and b, eta of beta
transform.df.parms = function(parms.df, out.format = c('b','beta','eta'), include.zero = TRUE)
{
  # start with many checks
  out.format = match.arg(out.format)
  colnames(parms.df) = tolower(colnames(parms.df))
  
  in.format = intersect(colnames(parms.df), c('b','beta','eta'))
  
  if(length(in.format) == 0)
    stop('parameters must contain at least one of following columns: b, beta, eta')
  
  if(length(in.format)>1)
  {
    in.format = in.format[1]
    message(paste0("Using '",in.format,"' as input parameter"))
  }
  
  if(!all(c('item_id','item_score') %in% colnames(parms.df)))
    stop('parameters must contain the columns: item_id, item_score')
  
  if(any(parms.df$item_score%%1 > 0))
    stop("column 'item_score' must be integer valued")
  
  parms.df = parms.df %>% 
    mutate(item_id = as.character(.data$item_id), item_score = as.integer(.data$item_score)) %>%
    arrange(.data$item_id, .data$item_score)
  
  
  mm = parms.df %>% 
    group_by(.data$item_id) %>% 
    summarise(min_score = min(.data$item_score), max_score = max(.data$item_score)) %>%
    ungroup()
  
  in.zero = any(mm$min_score == 0)
  
  if(in.zero && any(mm$min_score) != 0)
    stop("Either all item paramters should include a zero score parameter or none should")
  
  if(any(mm$max_score == 0))
    stop('All items should contain at least one non-zero score parameter')
  
  if(any(mm$min_score<0))
    stop("Negative scores are not allowed")

  if(in.format == 'b' && any(parms.df$b < 0))
      stop("A 'b' parameter cannot be negative, perhaps you meant to include a 'beta' parameter?")
  
  fl = parms.df %>%
    mutate(rn = row_number()) %>%
    group_by(.data$item_id) %>% 
    summarize(first = min(.data$rn), last=max(.data$rn)) %>%
    ungroup()
  
  args = list(first = fl$first, last = fl$last, parms.df = parms.df, 
              out.zero = include.zero, in.zero = in.zero)
  do.call(get(paste0(in.format,'2',out.format)), args)
}




##########################################
#' Fit the extended nominal response model
#'
#' Fits an Extended NOminal Response Model (ENORM) using conditional maximum likelihood (CML)
#' or a Gibbs sampler for Bayesian estimation.
#'
#'
#' @param dataSrc Data source: a dexter project db handle or a data.frame with columns: person_id, item_id, item_score
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param fixed_params Optionally, a prms object from a previous analysis or 
#' a data.frame with columns: item_id, item_score (omitting 0 score category) and beta. To facilitate the user in entering parameter values, we assume the parameterisation used by OPLM; in short, beta's are thresholds between categories. At this moment, it is not possible to fix some but not all categories of an item.
#' @param method If CML, the estimation method will be Conditional Maximum Likelihood;
#' otherwise, a Gibbs sampler will be used to produce a sample from the posterior
#' @param nIterations Number of Gibbs samples when estimation method is Bayes. The maximum 
#' number of iterations when using CML.
#' @return An object of type \code{prms}. The prms object can be cast to a data.frame of item parameters 
#' using function `coef` or used directly as input for other Dexter functions.
#' 
#' @references 
#' Maris, G., Bechger, T.M. and San-Martin, E. (2015) A Gibbs sampler for the (extended) marginal Rasch model. 
#' Psychometrika. 2015; 80(4): 859–879. 
#' 
#' @seealso functions that accept a prms object as input: \code{\link{ability}}, \code{\link{plausible_values}}
#'
fit_enorm = function(dataSrc, predicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), nIterations=500)
{
  method = match.arg(method)
  qtpredicate = eval(substitute(quote(predicate)))
  check_arg(dataSrc, 'dataSrc')
  check_arg(nIterations, 'integer', .length=1)
  env=caller_env()
  fit_enorm_(dataSrc, qtpredicate = qtpredicate, fixed_params = fixed_params, method=method, nIterations=nIterations, env=env)
}


fit_enorm_ = function(dataSrc, qtpredicate = NULL, fixed_params = NULL, method=c("CML", "Bayes"), nIterations=500, env=NULL) 
{
  method = match.arg(method)
  if(is.null(env)) env = caller_env()
  
  r = get_resp_data(dataSrc, qtpredicate, summarised=FALSE, env=env)
  
  x = r$x
  design = r$design
  if(nrow(x) == 0) stop('no data to analyse')
  if(nrow(design) == 1) stop('There are responses to only one item in your selection, this cannot be calibrated.') 
  
  if(length(unique(design$booklet_id)) > 1)
  {
    im = as.matrix(table(design$item_id, design$booklet_id))
    wm = crossprod(im, im)
    diag(wm) = 0
    if(!design_is_connected(list(im=im, wm=wm))) stop('Your design is not connected')  
  }
  
  itm_max = x %>% 
    group_by(.data$item_id) %>% 
    summarise(maxScore = max(.data$item_score)) %>% #, nsc = n_distinct(.data$item_score)) %>%
	  ungroup()
  
  if(any(itm_max$maxScore == 0)) 
    stop('One or more items has a maximum score of 0')
  # for now, just error. May possibly become warning later
  

  ssBIS = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$item_score) %>% 
    summarise(sufI=n(), sufC=sum(.data$item_score * .data$sumScore)) %>% 
    ungroup()
  

  plt = x %>% 
    group_by(.data$booklet_id, .data$item_id, .data$sumScore) %>% 
    summarise(meanScore=mean(.data$item_score), N=n()) %>% 
    ungroup()
  
  # max booklet scores
  maxScores = itm_max %>%
    inner_join(design, by='item_id') %>%
    group_by(.data$booklet_id) %>%
    summarise(maxTotScore = sum(.data$maxScore))
  
  # booklets 0:maxscore
  allScores = maxScores %>% 
    group_by(.data$booklet_id) %>%
    do(tibble(sumScore=0:.$maxTotScore)) %>%
    ungroup()
  

  stb = distinct(x,.data$person_id,.data$booklet_id,.keep_all=T) %>%
    group_by(.data$booklet_id, .data$sumScore) %>%
    summarize(N=n()) %>%
    right_join(allScores, by=c('booklet_id','sumScore')) %>%
    mutate(N=coalesce(.data$N, 0L)) %>%
    arrange(.data$booklet_id, .data$sumScore)

  ssIS = ssBIS %>% 
    group_by(.data$item_id, .data$item_score) %>%
    summarise(sufI = sum(.data$sufI)) %>%
    ungroup() %>%
    arrange(.data$item_id, .data$item_score)
  
  ssI  = ssIS %>% 
    group_by(.data$item_id) %>%
    summarise(nCat = n()) %>% 
    mutate(first = cumsum(.data$nCat) - .data$nCat + 1L,last = cumsum(.data$nCat)) %>%
    ungroup() %>%
    arrange(.data$item_id)
  
  if(any(ssI$nCat == 1)) stop('One or more items are without score variation')
  
  design = design %>%
    inner_join(ssI,by='item_id') %>%
    arrange(.data$booklet_id, .data$first)
  
  a = ssIS$item_score
  
  bkl = lapply(split(design,design$booklet_id),
    function(bkd)
    {
        booklet_id = bkd$booklet_id[1]
        scoretab = stb$N[stb$booklet_id==booklet_id]
        list(booklet_id = booklet_id,
             first = bkd$first,
             last = bkd$last,
             scoretab = scoretab,
             m = sum(scoretab),
             lambda = rep(1,length(scoretab)))
    })

  
  #it_sc_lab = paste0(ssIS$item_id[-ssI$first], "_",ssIS$item_score[-ssI$first])
  fixed_b = NULL
  has_fixed_parms = !is.null(fixed_params)
  
  ## deal with fixed parameters
  if(has_fixed_parms)
  {
    if(inherits(fixed_params,'prms'))
    {
      if (fixed_params$inputs$method!="CML")
        message("Posterior means are taken as values of fixed parameters")
      
      fixed_params = fixed_params$inputs$ssIS %>%
        add_column(b = if.else(fixed_params$inputs$method=="CML", fixed_params$est$b, colMeans(fixed_params$est$b))) %>%
        right_join(ssIS, by=c('item_id','item_score')) %>%
        arrange(.data$item_id,.data$item_score) 
    } else
    {
      # transform the fixed params to the b parametrization dexter uses internally
      fixed_params = transform.df.parms(fixed_params, out.format = 'b', include.zero = TRUE) 
    }  
    # check for missing categories in fixed_params, necessary?
    missing_cat = ssIS %>% 
      semi_join(fixed_params, by='item_id') %>%
      left_join(fixed_params, by=c('item_id','item_score')) %>%
      filter(is.na(.data$b) & .data$item_score != 0) 
      
    if(nrow(missing_cat) > 0)
    {
      cat(paste('Some score categories are fixed while some are not, for the same item.',
                'Dexter does not know how to deal with that.\nThe following score categories are missing:\n'))
      missing_cat %>% 
        select(.data$item_id, .data$item_score) %>%
        arrange(.data$item_id, .data$item_score) %>%
        as.data.frame() %>%
        print()
      stop('missing score categories for fixed items, see output')
    }
      
    fixed_b = fixed_params %>%
      right_join(ssIS, by=c('item_id','item_score')) %>%
      arrange(.data$item_id,.data$item_score) %>%
      pull(.data$b)
    
    # this test fails for some reason
    if(!anyNA(fixed_b)) stop('nothing to calibrate, all parameters are fixed')
  }
  
  if (method=="CML"){
    result = try(calibrate_CML(booklet=bkl, sufI=ssIS$sufI, a=ssIS$item_score, 
                               first=ssI$first, last=ssI$last, nIter=nIterations,
                               fixed_b=fixed_b))
    #names(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score)
    #row.names(result$beta)=it_sc_lab
    #rownames(result$acov.beta)=it_sc_lab
    #colnames(result$acov.beta)=it_sc_lab
  } else 
  {
    design = design %>%
      mutate(bn=dense_rank(.data$booklet_id)) %>%
      arrange(.data$first, .data$bn)
    
    itemListInt = split(design$bn, design$item_id)
    
    b = exp(runif(nrow(ssIS), -1, 1))
    
    result = try(calibrate_Bayes(itemList=itemListInt, booklet=bkl, sufI=ssIS$sufI, b=b, a=a, 
                                 first=ssI$first, last=ssI$last, nIter=nIterations,
                                 fixed_b=fixed_b))
    #colnames(result$b)= paste0(ssIS$item_id, "_",ssIS$item_score) 
    #colnames(result$beta)=it_sc_lab
  }
  if (inherits(result, "try-error")) stop('Calibration failed')
  
  abl_tables = new.env(parent = emptyenv())
  abl_tables$mle = NULL
  
  outpt = list(est=result, 
               inputs=list(bkList=bkl, ssIS=ssIS, ssI=ssI, design=r$design,plt=plt,
                            method=method, has_fixed_parms = has_fixed_parms), 
               abl_tables = abl_tables,
               xpr=deparse(qtpredicate))
  class(outpt) = append('prms', class(outpt)) 
  outpt
}




#' Plot for the extended nominal Response model
#' 
#' The plot shows 'fit' by comparing the expected score based on the model (grey line)
#' with the average scores based on the data (black line with dots) for groups of students
#' with similar estimated ability.
#' 
#' @param x object produced by fit_enorm
#' @param item_id which item to plot, if NULL, one plot for each item is made
#' @param nbins number of ability groups
#' @param ci confidence interval for the error bars, between 0 and 1. 0 means no error bars.
#' Default = 0.95 for a 95\% confidence interval
#' @param ... further arguments to plot
#' 
#' @method plot prms
#' 
plot.prms = function(x,item_id=NULL, nbins=5, ci = .95, ...)
{
  if(is.null(x$inputs$plt))
  {
    if(inherits(x,'mst_enorm')) stop('Sorry, the plot method for enorm_mst will only be supported in dexterMST from version 0.1.1 onwards')
    stop('Sorry, the plot method is only available for parameter objects produced with dexter 0.8.1 or later')
  }
  
  # if no item id provided plot them all
  if(is.null(item_id))
    item_id = x$inputs$ssI$item_id

  if(is.null(x$abl_tables$mle))
    x$abl_tables$mle = ability_tables(x, standard_errors=FALSE, method='MLE')
  
  if(length(item_id) > 1)
  {
    for(item in item_id) plot(x, item_id=item, nbins=nbins, ci=ci, ...)
    
    return(invisible(NULL))
  }
  # for dplyr
  item_id_ = item_id
  
  expf = expected_score(x, items = item_id)

  max_score = x$inputs$ssIS %>%
    filter(.data$item_id == item_id_) %>%
    pull(.data$item_score) %>%
    max()
  
  plt = x$inputs$plt %>%
    filter(.data$item_id==item_id_) %>%
    inner_join(x$abl_tables$mle, by=c('booklet_id','sumScore')) %>%
    filter(is.finite(.data$theta)) %>%
    mutate(abgroup = weighted_ntile(.data$theta, .data$N, n = nbins)) %>%
    group_by(.data$abgroup) %>%
    summarize(gr_theta = weighted.mean(.data$theta,.data$N), avg_score = weighted.mean(.data$meanScore,.data$N), n=sum(.data$N)) %>%
    ungroup() %>%
    mutate(expected_score = expf(.data$gr_theta))
  
  rng = max(plt$gr_theta) - min(plt$gr_theta)
  rng = c(min(plt$gr_theta)-.5*rng/nbins,
          max(plt$gr_theta)+.5*rng/nbins)
  
  plot.args = merge_arglists(list(...),
                             default=list(bty='l',xlab = expression(theta), ylab='score',main=item_id),
                             override=list(x = rng,y = c(0,max_score), type="n"))
  
  plot.args$main = fstr(plot.args$main, list(item_id=item_id))
  plot.args$sub = fstr(plot.args$sub, list(item_id=item_id))
  
  do.call(plot, plot.args)
  lines(plt$gr_theta,plt$expected_score, col='grey80') 
  
  if(!is.null(ci) && !is.na(ci) && ci !=0)
  {
    if(ci>1)
    {
      if(ci<100) ci = ci/100
    }  
    if(ci<0 || ci >= 1)
      stop('confidence interval must be between 0 and 1')
    
    qnt = abs(qnorm((1-ci)/2))
    
    cmin = function(p, n) pmax(0, p - qnt * sqrt(p*(1-p)/n))
    cmax = function(p, n) pmin(1, p + qnt * sqrt(p*(1-p)/n))
    
    arrows(plt$gr_theta, max_score*cmin(plt$expected_score/max_score, plt$n), 
           plt$gr_theta, max_score*cmax(plt$expected_score/max_score, plt$n), 
           length=0.05, angle=90, code=3, col='grey80')
  }
  
  lines(plt$gr_theta,plt$avg_score)  
  points(plt$gr_theta, plt$avg_score)
  invisible(NULL)
}



print.prms = function(x, ...){
  p = paste0( 'Parameters for the Extended Nominal Response Model\n\n',
              'Method: ', x$inputs$method, ', ',
              ifelse(x$inputs$method == 'CML',
                     paste0('converged in ',x$est$n_iter, ' iterations'),
                     paste0('number of Gibbs samples: ',nrow(x$est$beta))),
              '\nitems: ', nrow(x$inputs$ssI), 
              '\nresponses: ', sum(x$inputs$ssIS$sufI),'\n\n',
              'Use coef() or coefficients() to extract the item parameters.\n')
    
  cat(p)
  invisible(x)
}



#' extract enorm item parameters
#' 
#' @param object an enorm parameters object, generated by the function \code{\link{fit_enorm}}
#' @param bayes_hpd_b width of Bayesian highest posterior density interval around mean_beta, 
#'  value must be between 0 and 1, default is 0.95 
#' @param ... further arguments to coef are ignored
#'  
#' @return 
#' Depends on the calibration method:
#' \describe{
#' \item{for CML}{a data.frame with columns: item_id, item_score, beta, SE_beta}
#' \item{for Bayes}{a data.frame with columns: item_id, item_score, mean_beta, SD_beta, -bayes_hpd_b_left-, -bayes_hpd_b_right-}
#' }
#' 
#' 
#' 
coef.prms = function(object, bayes_hpd_b = 0.95, ...)
{

  
  if(bayes_hpd_b <= 0 ||  bayes_hpd_b >= 1)
    stop('args$bayes_hpd_b must be between 0 and 1')
  
  x = object
  
  hpd=function(x, conf=args$bayes_hpd_b){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return(data.frame(l=x[ nnn ],r=x[ n-nn+nnn ]))
  }
  
  if (x$inputs$method=="CML")
  {
    atab=data.frame(item_id=x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    item_score=as.integer(x$inputs$ssIS$item_score[-x$inputs$ssI$first]),
                    beta=x$est$beta,
                    SE_beta=sqrt(diag(x$est$acov.beta)),stringsAsFactors=FALSE)
  }
  
  if (x$inputs$method=="Bayes"){
    hh = map_df(apply(x$est$beta,2,hpd),c)
    atab=data.frame(item_id = x$inputs$ssIS$item_id[-x$inputs$ssI$first],
                    a = x$inputs$ssIS$item_score[-x$inputs$ssI$first],
                    mb = colMeans(x$est$beta),
                    sdb = apply(x$est$beta, 2, sd),
                    hpdl = hh[,1], hpdr=hh[,2],stringsAsFactors=FALSE)
    colnames(atab)=c("item_id" ,"item_score", "mean_beta", "SD_beta", 
                     sprintf("%i_hpd_b_left", round(100 * bayes_hpd_b)),
                     sprintf("%i_hpd_b_right", round(100 * bayes_hpd_b)))
  }
  rownames(atab) = NULL
  return(atab)
}


#' Functions of theta
#' 
#' returns information function or expected score function for a single item, an arbitrary group of items or all items
#' 
#' @param parms object produced by fit_enorm
#' @param items vector of one or more item_id's. If NULL and booklet_id is also NULL, all items in parms are used
#' @param booklet_id id of a single booklet (e.g. the test information function), if items is not NULL this is ignored
#' @param which.draw the number of the random draw (only applicable if calibration method was Bayes). If NULL, the mean 
#' beta parameter will be used
#' 
#' @return a function which accepts a vector of theta's and returns an equal length vector 
#' with the information estimate or the expected score
#' 
#' @examples
#' 
#' db = start_new_project(verbAggrRules,':memory:')
#' add_booklet(db,verbAggrData, "agg")
#' p = fit_enorm(db)
#' 
#' # plot information function for single item
#' 
#' ifun = information(p, "S1DoScold")
#' 
#' plot(ifun,from=-4,to=4)
#' 
#' # compare test information function to the population ability distribution
#' 
#' ifun = information(p, booklet="agg")
#' 
#' pv = plausible_values(db,p)
#' 
#' op = par(no.readonly=TRUE)
#' par(mar = c(5,4,2,4))
#' 
#' plot(ifun,from=-4,to=4, xlab='theta', ylab='test information')
#' 
#' par(new=TRUE)
#' 
#' plot(density(pv$PV1), col='green', axes=FALSE,xlab=NA, ylab=NA,main=NA)
#' axis(side=4)
#' mtext(side = 4, line = 2.5, 'population density (green)')
#' 
#' par(op)
#' close_project(db)
#' 
information = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='information')
}

#' @rdname information
expected_score = function(parms, items=NULL, booklet_id=NULL, which.draw=NULL)
{
  theta_function(parms, items=items, booklet=booklet_id, which.draw=which.draw, what='expected')
}
# TO DO: should be able to accept a data.frame as wel as parms in the future
# escore can be more efficient with the main loop in C
theta_function = function(parms, items=NULL, booklet=NULL, which.draw=NULL, what=c('information','expected'))
{
  what = match.arg(what)
  if(is.null(items) && is.null(booklet))
  {
    fl = parms$inputs$ssI
  } else if(!is.null(items))
  {
    items = unique(items)
    
    if(length(setdiff(items,parms$inputs$ssI$item_id))>0)
    {
      message('unknown items:')
      print(sort(setdiff(items,parms$inputs$ssI$item_id)))
      stop('Some items were not found, see output')
    }
    
    fl = parms$inputs$ssI %>% 
      semi_join(tibble(item_id=items),by='item_id') %>% 
      arrange(.data$first)
    
  } else
  {
    booklet = as.character(booklet)
    check_arg(booklet,'character',.length=1)
    
    fl = parms$inputs$design %>%
      filter(.data$booklet_id==booklet) %>%
      inner_join(parms$inputs$ssI, by='item_id') %>%
      arrange(.data$first)
    if(nrow(fl)==0)
      stop(paste0('booklet `',booklet,'` not found'))
  }
  
  if(is.matrix(parms$est$b))
  {
    if(is.null(which.draw))
    {
      parms$est$b = colMeans(parms$est$b)
    } else
    {
      if(which.draw<1 || which.draw > nrow(parms$est$b))
        stop('parameter `which.draw` out of range')
      parms$est$b = as.vector(parms$est$b[which.draw,])
    }
  }
  if(what=='information')
  {
    out = function(theta)
    {
      check_arg(theta,'numeric')
      if(any(is.na(theta) | is.nan(theta)))
        stop('theta may not contain nan/NA values')  
      
      res = rep(0,length(theta))
      res[is.finite(theta)] = 
        IJ_(parms$est$b,parms$inputs$ssIS$item_score,fl$first, fl$last, theta[is.finite(theta)])$I
      # extremely large values overflow to NaN, recover as 0
      res[is.nan(res)] = 0
      res
    }
    class(out) = append('inf_func',class(out))
    
  } else if(what=='expected')
  {
    max_score = sum(parms$inputs$ssIS$item_score[fl$last])
    out = function(theta)
    {
      check_arg(theta,'numeric')
      if(any(is.na(theta) | is.nan(theta))) 
        stop('theta may not contain nan/NA values') 
      
      res = rep(0,length(theta))
      
      res[is.finite(theta)] =
        sapply(theta[is.finite(theta)], E_score, b=parms$est$b, a=parms$inputs$ssIS$item_score, first=fl$first, last=fl$last)
      
      res[is.infinite(theta) & theta > 0] = max_score
      # extremely large values of theta overflow to NaN (small values undeflow to zero, which is fine)
      res[is.nan(res)] = max_score
      res
    }
    class(out) = append('exp_func',class(out))
  }
  out
}

print.inf_func = function(x,...) cat('Information function I(theta)\n')
print.exp_func = function(x,...) cat('Expected score function E(theta)\n')

