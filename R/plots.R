

# colors derived from http://colorbrewer2.org
qcolors = function(n, user_colors=NULL)
{
  if(is.function(user_colors))
  {
    pal = user_colors(n)
  } else if(length(user_colors) > 1 || length(user_colors) == n)
  {
    pal = user_colors
  } else
  {
    pal = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999",
            "#BC80BD", "#CCEBC5","#FFED6F","#FB8072","#80B1D3","#FDB462","#8DD3C7","#FFFFB3","#BEBADA",
            "#B3DE69","#A50F15","#08306B","#00441B","#54278F","#fc4E2A","#525252","#66C2A4")
  }
  rep_len(pal,n)
}

dvcolors = function(n, user_colors=NULL)
{
  if(is.function(user_colors))
  {
    pal = user_colors(n)
  } else if(length(user_colors) > 1 || length(user_colors) == n)
  {
    pal = user_colors
  } else
  {
    pal = c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5',
            '#3288BD','#5E4FA2')
  }
  colorRampPalette(pal)(n)
}



lighten = function(color, factor = 0.5) {
  if ((factor > 1) || (factor < 0)) stop("factor needs to be within [0,1]")
  col = col2rgb(color)
  col = col + (255 - col)*factor
  col = rgb(t(col), maxColorValue=255)
  col
}


draw_curtains = function(qnt)
{
  if(!is.null(qnt))
  {
    usr = par('usr')
    rect(usr[1], usr[3], qnt[1], usr[2], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
    rect(qnt[2], usr[3], usr[2], usr[4], col="#EEEEEE", border=NA,xpd=FALSE,lwd=0)
    # rect will sometimes cover the axis (redrawing the axis would probably solve this too)
    abline(h=usr[3])
    abline(v=usr[1])
  }
}


###################################################
#' Distractor plot
#'
#' Produce a diagnostic distractor plot for an item
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: person_id, item_id, response, item_score
#' and optionally booklet_id
#' @param item_id The ID of the item to plot. A separate plot will be produced
#' for each booklet that contains the item, or an error message if the item_id
#' is not known. Each plot contains a non-parametric regression of each possible
#' response on the total score.
#' @param predicate An optional expression to subset data, if NULL all data is used
#' @param legend logical, whether to include the legend. default is TRUE
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param adjust factor to adjust the smoothing bandwidth respective to the default value
#' @param col vector of colors to use for plotting
#' @param ... further arguments to plot.
#' @return 
#' Silently, a data.frame of response categories and colors used. Potentially useful if you want to customize the legend or 
#' print it separately
#' @details 
#' Customization of title and subtitle can be done by using the arguments main and sub. 
#' These arguments can contain references to the variables item_id, booklet_id, item_position(if available),
#' pvalue, rit and rir. References are made by prefixing these variables with a dollar sign. Variable names may be postfixed 
#' with a sprintf style format string, e.g. 
#' \code{distractor_plot(db, main='item: $item_id', sub='Item rest correlation: $rir:.2f')}
#' 
distractor_plot = function(dataSrc, item_id, predicate=NULL, legend=TRUE, curtains=10, adjust=1, col=NULL, ...){  
  check_dataSrc(dataSrc)

  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  item_id = as.character(item_id)
  check_string(item_id)
  item = item_id
  
  iprop = list()
  if(is_db(dataSrc))
    iprop = as.list(dbGetQuery_param(dataSrc,'SELECT * FROM dxItems WHERE item_id= :item;', 
                                     tibble(item=item)))
  
  
  if(is.null(qtpredicate) && is_db(dataSrc))
  {
    # pre process a little to make things faster
  	booklets = dbGetQuery_param(dataSrc,
  	     'SELECT booklet_id FROM dxbooklet_design WHERE item_id=:item;', tibble(item=item)) %>%
  		pull('booklet_id') %>%
  	  sql_quote("'") %>%
  		paste(collapse=',')
  	
    qtpredicate = sql(paste0("booklet_id IN(",booklets,")"), 'booklet_id')
  } 
  
  respData = get_resp_data(dataSrc, qtpredicate = qtpredicate, extra_columns='response', env=env, summarised=FALSE) %>%
    filter(.data$item_id == item, .recompute_sumscores = FALSE )
  

  if (nrow(respData$design) == 0) 
    stop(paste("Item", item, "not found in dataSrc."))
  

  if('item_position' %in% colnames(respData$design))
  {
    ipos = respData$design
  } else if(is_bkl_safe(dataSrc, qtpredicate, env) && is_db(dataSrc))
  {
    ipos = dbGetQuery(dataSrc, paste("SELECT booklet_id, item_position FROM dxbooklet_design WHERE item_id=", 
                                     sql_quote(item,"'")))
  } else if(inherits(dataSrc,'data.frame') && 'item_position' %in% colnames(dataSrc))
  {
    ipos = respData$design = dataSrc %>%
      filter(.data$item_id==item) %>%
      distinct(.data$booklet_id,.keep_all=TRUE)
  } else
  {
    ipos = respData$design
  }
  

  default.args = list(sub = "Pval: $pvalue:.2f, Rit: $rit:.3f, Rir: $rir:.3f", 
                      xlab = "Sum score", ylab = "Proportion", cex.sub = 0.8, xaxs="i", bty="l")
  
  default.args$main = ifelse('item_position' %in% colnames(ipos),
                             '$item_id, pos. $item_position in booklet $booklet_id',
                             '$item_id in booklet $booklet_id')
  
  user.args = list(...)
  
  rsp_counts = respData$x %>%
    count(.data$booklet_id, .data$response, .data$item_score, .data$booklet_score)

  max_score = max(rsp_counts$item_score)

  rsp_colors = rsp_counts %>%
    group_by(.data$response) %>%
    summarise(n = sum(.data$n)) %>%
    ungroup() %>%
    arrange(desc(.data$n), .data$response) 
  
  # cannot make this a list since a response can be an empty string
  rsp_colors = rsp_colors %>% add_column(color = qcolors(nrow(rsp_colors),col))
  
  
  lapply(split(rsp_counts, rsp_counts$booklet_id), function(y)
  {
    booklet = y$booklet_id[1]
    
    bkl_scores = y %>% 
      group_by(.data$booklet_score) %>% 
      summarise(n = sum(.data$n)) %>%
      ungroup() %>%
      arrange(.data$booklet_score)
    
    if(nrow(bkl_scores)>1)
    {
      N = sum(bkl_scores$n)
      
      labs = modifyList(iprop,
        list(pvalue = weighted.mean(y$item_score, y$n)/max_score,
                  rit = weighted_cor(y$item_score, y$booklet_score, y$n),
                  rir = weighted_cor(y$item_score, y$booklet_score - y$item_score, y$n),
                  n = N,
                  item_position = filter(respData$design, .data$booklet_id==booklet)$item_position,
                  booklet_id=booklet,
                  item_id=item))
      
      plot.args = merge_arglists(user.args, 
                                 default=default.args,
                                 override=list(x = c(0,max(y$booklet_score)), y = c(0,1), type="n"))
      
      plot.args$main = fstr(plot.args$main, labs)
      plot.args$sub = fstr(plot.args$sub, labs)
      
      qua = curtains/200
      qnt = NULL
      if(qua>0 && qua<.5) {
        #qnt = quantile(rep(bkl_scores$booklet_score, bkl_scores$n), c(qua,1-qua))
        qnt = weighted_quantile(bkl_scores$booklet_score, bkl_scores$n, c(qua,1-qua))
      }
      
      do.call(plot, plot.args)
      draw_curtains(qnt)
      
      dAll = density(bkl_scores$booklet_score, n = 512, weights = bkl_scores$n/N, adjust=adjust)
      
      lgnd = y %>% 
        group_by(.data$response)  %>% 
        do({
          k = rsp_colors[rsp_colors$response == .$response[1],]$color
          
          if(nrow(.)==1)
          {
            yval = .$n/sum(filter(y, .data$booklet_score==.$booklet_score)$n)
            points(.$booklet_score,yval,col=k,pch=16,xpd=TRUE)
          } else
          {
            dxi = density(.$booklet_score, weights = .$n/sum(.$n), n = 512,
                        bw = dAll$bw, from = min(dAll$x), to = max(dAll$x))
            yy = dxi$y/dAll$y * sum(.$n)/N
            lines(dAll$x, yy, col = k, lw = 2)
          }
          tibble(col = k, label = paste0(.$response[1]," (", .$item_score[1], ")"))
        })
      
      if(legend && NROW(lgnd)>0)
      {
        graphics::legend("topleft", legend = lgnd$label, inset=c(0.01, 0),
                         lty = 1, col = lgnd$col, cex = 0.8,lwd=2, box.lty = 0)
      }
    }
  })

  invisible(df_format(rsp_colors))
}

#' Profile plot
#'
#'
#' @param dataSrc a connection to a dexter database or a data.frame with columns: 
#' person_id, item_id, item_score and the item_property and the covariate of interest.
#' @param item_property The name of the item property defining the domains. 
#' The item property should have exactly two distinct values in your data
#' @param covariate name of the person property used to create the groups. 
#' There will be one line for each distinct value.
#' @param predicate An optional expression to filter data, if NULL all data is used
#' @param model "IM" (default) or "RM" where "IM" is the interaction model and 
#' "RM" the Rasch model. The interaction model is the default as it fits 
#' the data better or at least as good as the Rasch model.
#' @param x Which value of the item_property to draw on the x axis, if NULL, one is chosen automatically
#' @param col vector of colors to use for plotting
#' @param ... further arguments to plot
#' @details 
#' Profile plots can be used to investigate whether two (or more) groups of respondents 
#' attain the same test score in the same way. The user must provide a  
#' (meaningful) classification of the items in two non-overlapping subsets such that 
#' the test score is the sum of the scores on the subsets. 
#' The plot shows the probabilities to obtain 
#' any combinations of subset scores with thin gray lines indicating the combinations 
#' that give the same test score. The thick lines connect the most likely 
#' combination for each test score in each group.
#' When applied to educational test data, the plots can be used to detect differences in the 
#' relative difficulty of (sets of) items for respondents that belong to different 
#' groups and are matched on the test score. This provides a content-driven way to 
#' investigate differential item functioning. 
#'
#' @examples

#' db = start_new_project(verbAggrRules, ":memory:", 
#'                          person_properties=list(gender="unknown"))
#' add_booklet(db, verbAggrData, "agg")
#' add_item_properties(db, verbAggrProperties)
#' profile_plot(db, item_property='mode', covariate='gender')
#' 
#' close_project(db)
#' 
#' 
profile_plot = function(dataSrc, item_property, covariate, predicate = NULL, model = c("IM","RM"), x = NULL, col = NULL, ...) 
{
  check_dataSrc(dataSrc)
  check_string(item_property)
  check_string(covariate)
  model = match.arg(model)
  user.args = list(...)
  if(is_db(dataSrc))
  {
    item_property = dbValid_colnames(item_property)
    covariate = dbValid_colnames(covariate)
  }

  qtpredicate = eval(substitute(quote(predicate)))
  env = caller_env()
  
  if(inherits(dataSrc,'matrix'))
	stop('profile_plot does not accept a matrix dataSrc')
  
  respData = get_resp_data(dataSrc, qtpredicate, extra_columns = covariate, env = env)  %>%
	  intersection()
  
  respData$x[[covariate]] = ffactor(respData$x[[covariate]])
  
    
  if(is_db(dataSrc) && item_property %in% dbListFields(dataSrc,'dxitems'))
  {
    domains = dbGetQuery(dataSrc, paste("SELECT item_id,", item_property, "FROM dxitems;")) %>%
      semi_join(tibble(item_id = levels(respData$x$item_id)), by='item_id')
    
    domains$item_id = ffactor(domains$item_id, levels = levels(respData$x$item_id))
    
  } else
  {
    domains = distinct(dataSrc, .data$item_id, .data[[item_property]])
    
    if(nrow(domains) > n_distinct(domains$item_id))
      stop("some items belong to multiple domains, this is not allowed")
  }
  
  domains = arrange(domains, .data$item_id)
  
  props = unique(domains[[item_property]])
  if(length(props) != 2)
    stop('this function needs an item_property with 2 unique values in your data')
  
  # order props to get the x and y axis right
  if(!is.null(x) && x %in% props) 
    props = c(x, props[props!=x])
  

  AB = list(A = which(domains[[item_property]] == props[1]), 
            B = which(domains[[item_property]] == props[2]))
  
  # fit interaction model and compute ssTable
  models = by(respData, covariate, fit_inter_, regs=FALSE) 
  
  tt = lapply(models, function(m)
  {
      if(model=="IM")
      {
        SSTable(b=m$est$bIM, a=m$inputs$ssIS$item_score, 
                first=m$inputs$ssI$first, last=m$inputs$ssI$last, AB = AB, m$est$cIM)
      } else
      {
        SSTable(b=m$est$bRM, a=m$inputs$ssIS$item_score, 
                first=m$inputs$ssI$first, last=m$inputs$ssI$last, AB = AB)
      }
  })

  psbl = Reduce(union,lapply(models,function(m) m$est$possible_scores))
  
  maxA = max(sapply(tt, nrow)) - 1L
  maxB = max(sapply(tt, ncol)) - 1L

  default.args = list(main="Profile plot", xlab=props[1], 
                      ylab=props[2],xlim=c(0,maxA),ylim=c(0,maxB),bty='l')
  do.call(plot, 
          merge_arglists(user.args, 
                         default=default.args,
                         override=list(x=c(0,maxA), y=c(0,maxB),type="n")))
  clip(0,maxA,0,maxB) 
  segments(0L, psbl, psbl, 0L, col='gray', xpd=FALSE)
  text(pmin(maxA,psbl)+.4, pmax(0, psbl-maxA)-.4,psbl,cex=.6,col="lightgray", xpd=TRUE)
  
  colors = qcolors(length(tt), col)
  
  for (i in seq_along(tt)) 
  {
    lns = sapply(psbl, function(s){
        indx = cbind(1L+(s:0),1L+(0:s))
        indx = indx[indx[,1]<=nrow(tt[[i]]) & indx[,2] <= ncol(tt[[i]]),, drop=FALSE]
        indx = indx[which.max(tt[[i]][indx]),, drop=FALSE]
        if(nrow(indx) > 0 && tt[[i]][indx] > 0)
          return(indx)
        rep(NA_integer_, 2)
      }) %>%
      t() %>%
      apply(2, '-', 1) 
    
    lns = lns[!(is.na(lns[,1]) | is.na(lns[,2])),]
    lines(lns,col=colors[i], lw=2)
  }

  
  legend("topleft", legend=names(tt),lty=1, col=colors, cex=.7, 
         box.lty=0, bg='white')
  invisible(NULL)
}


