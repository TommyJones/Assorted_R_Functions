require(forecast)
require(tseries)
require(sandwich)
require(lmtest)
require(car)

Lag <- function(x, p){
  # inputs are x, a vector, and p, the number of periods to lag
  # x must be sorted by its time index beforehand and NAs should be inserted for any missing periods of data beforehand
  
  lag.x <- c(rep(NA, p), x[1:(length(x) - p)])
  names(lag.x) <- names(x) # if you are using names to index your time series, this ensures you have the right names
  
  return(lag.x)
  
}

MakeTs<-function(var, time.index){
  #these steps ensure that there is an observation (even if it is NA) for each time period
  time.index <- as.numeric(time.index)
  
  new.index <- data.frame(index=min(time.index, na.rm=TRUE):max(time.index, na.rm=TRUE))
  
  ts.data <- data.frame(index=time.index, var=var)
  
  ts.data <- merge(ts.data, new.index, ALL=TRUE)
  
  #creat the time series object
  ts.out <- ts(data=ts.data$var, start=min(ts.data$index, na.rm=TRUE), end=max(ts.data$index, na.rm=TRUE))
  
  return(ts.out)
}

GetDiff <- function(x, pct=FALSE, laplace.eps=0){
  # calculates period on period changes in x
  # if pct = FALSE (default) then level differences are returned
  # if pct = TRUE then percent changes are returned
  # eps is value to be used for laplace smoothing, if desired
  # x must be sorted by its time index beforehand and NAs should be inserted for any missing periods of data beforehand
  
  d <- x[ 2:length(x) ] - x[ 1:(length(x) - 1) ]
  result <- c(NA, d)
  names(result) <- names(x) # if you are using names to index your time series, this ensures you have the right names
  
  if( ! pct ){ return( result ) }
  
  else{
    p <- d / (x[ 1:(length(x) - 1) ] + laplace.eps) * 100 # uses laplace smoothing to not get infinite values
    result <- c(NA, p)
    names(result) <- names(x) # if you are using names to index your time series, this ensures you have the right names
    return( result )
  }
}


ChowTest <- function(y, x, restriction){
  # performs a chow test: h0: no break over obeservations denoted by restriction
  # y is the left hand side variable
  # x is a matrix of the right hand side variables
  # restriction is a logical denoting observations on either side of a break
  
  r <- lm(y ~ x)
  u1 <- lm(y[restriction] ~ x[restriction, ])
  u2 <- lm(y[!restriction] ~ x[!restriction, ])
  
  ssr.r <- sum( r$resid * r$resid )
  ssr.u <- sum( (u1$resid * u1$resid) ) + sum( (u2$resid * u2$resid) )
  
  k <- ncol(x) # number of regressors
  n <- nrow(x)
  df <- n - 2 * k
  
  chow <- ( ( ssr.r - ssr.u ) / k ) / (ssr.u / df )
  
  p.val <- 1-pf(chow, k, df)
  
  result <- list(chow=chow, p.val=p.val)
  
  return(result)
  
}



DummySat <- function(tseries, search.periods=names(tseries), incl.vars=NULL, alpha, data, method="forward"){
  # Takes a vector ordered by time period, performs two step dummy saturation, and returns a logical vector for those periods deemed significant
  # tseries is a string denoting the vector to be tested from data, its names must be indexed by the time period
  # search.periods is a vecotr of the same length as tseries giving the time index,
  # incl.vars is vector of strings denoting the names of variables to include for every regression from data. 
  # alpha is the p-value deemed to be "significant"
  # data is a data frame or matrix of data for your analysis. Its rownames must be indexed by the time period.
  # method is either "forward" or "backward". It determines the direction of asymetric bucketing. 
  # "forward" saturates the first third of search.periods and then the back two thirds in its two step regressions.
  # "backward" saturates the first two thirds and then the back one third of search.periods.
  
  tseries <- data[ , tseries ]
  
  # set up period dummies to test
  period.mat <- sapply(search.periods, function(p){
    return( as.numeric( rownames(data) %in% p ) )
  })
  colnames(period.mat) <- paste("p.", search.periods, sep="")
  rownames(period.mat) <- colnames( rownames(data) )
  
  third.mark <- ceiling(length(search.periods) / 3)
  
  if( method == "backward" ){
    periods.1.a <- colnames(period.mat[ , 1:(2 * third.mark) ])
    periods.2.a <- colnames(period.mat[ , (2 * third.mark + 1):ncol(period.mat) ])
  }else{ 
    periods.1.a <- colnames(period.mat[ , 1:third.mark ])
    periods.2.a <- colnames(period.mat[ , (third.mark + 1):ncol(period.mat) ])
  }
  # set up matrices of right hand side variables
  if( ! is.null(incl.vars) ){ 
    x.1.a <- cbind(data[ , incl.vars ], period.mat[ , periods.1.a ])
    x.2.a <- cbind(data[ , incl.vars ], period.mat[ , periods.2.a ])
    
  }else{ 
    x.1.a <- period.mat[ , periods.1.a ]
    x.2.a <- period.mat[ , periods.2.a ]
    
  }
  
  # perform initial regressions for part a
  reg.1.a <- lm(tseries ~ x.1.a)
  sum.reg.1.a <- summary(reg.1.a)
  
  reg.2.a <- lm(tseries ~ x.2.a)
  sum.reg.2.a <- summary(reg.2.a)
  
  # rename coefficient names in regression summary for easy identification
  rownames(sum.reg.1.a$coeff) <- c("Intercept", colnames(x.1.a))
  rownames(sum.reg.2.a$coeff) <- c("Intercept", colnames(x.2.a))
  
  # store significant impulse variables
  keep.1.a <- sum.reg.1.a$coeff[rownames(sum.reg.1.a$coeff) %in% periods.1.a , 4 ] <= alpha
  keep.2.a <- sum.reg.2.a$coeff[rownames(sum.reg.2.a$coeff) %in% periods.2.a , 4 ] <= alpha
  
  # perform our third regression for part a
  if( ! is.null(incl.vars) ){ 
    x.3.a <- cbind(incl.vars, period.mat[ , c(periods.1.a[ keep.1.a ], periods.2.a[ keep.2.a ]) ])
    colnames(x.3.a) <- c(colnames(incl.vars), c(periods.1.a[ keep.1.a ], periods.2.a[ keep.2.a ]))
  }else{ 
    x.3.a <- as.matrix( period.mat[ , c(periods.1.a[ keep.1.a ], periods.2.a[ keep.2.a ]) ] )
    colnames(x.3.a) <- c(periods.1.a[ keep.1.a ], periods.2.a[ keep.2.a ])
  }
  
  reg.3.a <- lm(tseries ~ x.3.a )
  sum.reg.3.a <- summary(reg.3.a)
  
  rownames(sum.reg.3.a$coeff) <- c("Intercept", colnames(x.3.a))
  
  
  # store final significant impulse variables
  sig.periods.a <- rownames(sum.reg.3.a$coeff)[rownames(sum.reg.3.a$coeff) %in% colnames(period.mat) & sum.reg.3.a$coeff[ , 4 ] <= alpha ]
  
  result.a <- search.periods[ paste("p.", search.periods, sep="") %in% sig.periods.a]
  
  return(result.a)
  
}

MakeStep <- function(steps, name){
  # steps is a list of step periods to try in start.year:2012
  # name the cohort or other identifier
  # returns a data frame indexed by start.year to 2012
  # column names of data frame are given by "names.YY.YY" where YY and YY give last two digets of start and end year of step
  
  result <- data.frame(sapply(steps, function(S){
    step <- as.numeric(start.year:2012 %in% S)
  }))
  
  names(result) <- sapply(steps, function(S){
    paste(name, min(S), max(S), sep=".")
  })
  
  result$Year <- start.year:2012
  
  result <- result[ , c("Year", names(result)[ names(result) != "Year" ]) ]
  
  return(result)
  
}

#########################################
# Make a function for getting differences
# out of panel data
#########################################
PanelDiff <- function(varnames, dframe, pct=FALSE){
  # returns a data frame of differences indexed by year and state
  tmp <- dframe[ , c("Year", "StateCode", varnames) ]
  
  tmp <- tmp[ order(tmp$StateCode, tmp$Year) , ]
  
  difflist <- lapply(unique(tmp$StateCode), function(state){
    statesub <- tmp[ tmp$StateCode == state , ]
    
    all.yrs <- min(statesub$Year,na.rm=TRUE):max(statesub$Year, na.rm=TRUE)
    all.yrs <- data.frame(all.yrs = all.yrs, stringsAsFactors=FALSE)
    
    statesub <- merge(statesub, all.yrs, by.x="Year", by.y="all.yrs")
    
    statesub[ , varnames ] <- apply(statesub[ , varnames ], 2, function(x) GetDiff(x, pct=pct) )
    
    return(statesub)
  })
  
  result <- do.call(rbind, difflist)
  return(result)
  
}


NwReg <- function(model.object, data){
  # updates regression based on particular data set (can be used for recursive regs)
  # calculates Newey West standard errors (1987)
  # returns a list with a regression table, newey west variance-covariance matrix, 
  # and lm object of new fit
  
  fit <- update(model.object, data=data, na.action=na.omit)
  
  nw.cov <- NeweyWest(fit, prewhite=0, lag=1) # newey west variance covariance matrix of coefficients
  
  ctest <- coeftest(x=fit, vcov.=nw.cov) # summary regression with nw std errors
  
  result <- data.frame(ctest[ , 1:4])
  names(result) <-  colnames(ctest)
  rownames(result) <- rownames(ctest)
  
  result <- list(result, nw.cov, fit)
  names(result) <- c("table", "nw.cov", "lm.object")
  
  return(result)
}

GetNwCi <- function(nwfit, conf.lvl){
  # takes result of NwReg() (nwfit) and returns a conf.lvl*100 % confidence interval for the coefficients of nwfit
  
  result <- data.frame(coeff=nwfit$table$Estimate)
  alpha <- 1 - conf.lvl
  
  ci <- nwfit$table$"Std. Error"*qt(p=1 - alpha/2, 
                                    (df=length(nwfit$lm.object$resid) - nrow(nwfit$table) ))
  
  result$high <- nwfit$table$Estimate + ci
  result$low <- nwfit$table$Estimate - ci
  
  return(result)
}


NwRecursive <- function(p.range, fit, dframe, time.index.name, conf.lvl){
  # recursive regressions with NW standard errors and confidence level
  # Uses NwReg() and NwCi()
  
  tmp <- lapply(p.range, function(end.year){
    data <- dframe[ dframe[ , time.index.name] <= end.year , ]
    
    try(nwfit <- NwReg(fit, data))
    
    if(exists("nwfit")){
      nw.ci <- GetNwCi(nwfit=nwfit, conf.lvl=conf.lvl)
    }else{
      nw.ci <- data.frame(coeff=rep(NA,length(fit$coefficients)), 
                          high=rep(NA,length(fit$coefficients)), 
                          low=rep(NA,length(fit$coefficients))
                          )
    }

    return(nw.ci)
    
  })
  
  names(tmp) <- p.range
  
  result <- lapply(1:length(fit$coefficients), 
                   function(x){
                     do.call(rbind, lapply(tmp, function(year) year[ x , ] ))
                   }
  )
  names(result) <- names(fit$coefficients)
  return(result)
}  

PanelApply <- function(x.frame, group.var, time.var, FUN){
  # Appies time series functions to panel data
  # input is a data frame of variables to which you want to apply FUN, indexed by time and group
  # output is a data frame with FUN applied indexed by time and group
  # FUN can be user defined or an out of the box function
  
  x.frame <- x.frame[ order(x.frame[ , group.var ], x.frame[ , time.var ]) , ]
  
  result <- lapply(unique(x.frame[ , group.var ]), 
                          function(GROUP){
                            
                            vars.by.group <- x.frame[ x.frame[ , group.var ] == GROUP , ]
                            
                            complete.periods <- min(vars.by.group[ , time.var ], na.rm=TRUE):max(vars.by.group[ , time.var ], na.rm=TRUE)
                            
                            complete.periods <- data.frame(complete.periods=complete.periods)
                            
                            vars.by.group <- merge(vars.by.group, complete.periods, by.x=time.var, by.y="complete.periods", all=TRUE)
                            
                            vars.by.group[ , setdiff(names(vars.by.group), c(time.var, group.var))] <- 
                              apply( data.frame(vars.by.group[ , setdiff(names(vars.by.group), c(time.var, group.var))]), 
                                     2, 
                                     FUN
                                     )
                            
                            return(vars.by.group[ ! is.na(vars.by.group[ , group.var ]) , ])
                          }
                   )
                   
    
  result <- do.call(rbind, result)  
    
  
  return(result) 
}
