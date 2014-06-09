require(randomForest)
require(compiler)

Rf <- cmpfun(function(data, class.var, num.predictors=floor(sqrt(ncol(data)-1)),  
                      num.models=min( 1000, choose( (ncol(data)-1), num.predictors) ), ...){
  result <- randomForest(x=data[ , setdiff(names(data), class.var) ], y=data[ , class.var ], mtry=num.predictors, ntrees=num.models, 
                         na.action=na.omit, replace=TRUE
                         )
  return(result)
})

RfEffectSize <- cmpfun(function( data, class.var, class.lvl, effect.var, mtry=4,
                                 ntree=min( 1000, choose( (ncol(data)-1), mtry)), ... ){
  result <- lapply( unique( data[ !is.na( data[ , effect.var]) , effect.var ]) , function(x){
#     print(x)
    class.data <- data[ data[ , effect.var ] == x , setdiff(names(data), effect.var) ]
    
    fit.1 <- randomForest(x=class.data[ , setdiff(names(class.data), class.var)], y=class.data[ , class.var], mtry=mtry ,ntree=ntree, replace=TRUE)
    
    class.data <- data[ data[ , effect.var ] != x , setdiff(names(class.data), effect.var) ]
    
    fit.2 <- randomForest(x=class.data[ , setdiff(names(class.data), class.var)], y=class.data[ , class.var], mtry=mtry ,ntree=ntree, replace=TRUE)
    
    pred.1 <- fit.1$votes[ , class.lvl ]
    pred.2 <- fit.2$votes[ , class.lvl ]
    
    diff <- mean(pred.1,na.rm=TRUE) - mean(pred.2, na.rm=TRUE)
    size <- diff
 
    mydist <- sapply( 1:5000, function(rep){
      p1 <- sample(pred.1, size=length(pred.1), replace=TRUE)
      p2 <- sample(pred.2, size=length(pred.1), replace=TRUE)
      
      result <- mean(p1, na.rm=TRUE) - mean(p2, na.rm=TRUE)
      
      return(result)
    })
    
    conf.95 <- quantile(mydist, c(0.025, 0.975), na.rm=TRUE)
    
    
    result <- data.frame(var=effect.var,
                         value=x,
                         size=size,
                         conf.95.low=conf.95[1],
                         conf.95.high=conf.95[2],
                         sig= (conf.95[1] > 0 & conf.95[2] > 0) | (conf.95[1] < 0 & conf.95[2] < 0) ,
                         stringsAsFactors=FALSE
                         )
    result[ , 1:2 ] <- apply(result[ , 1:2 ], 2, as.character)
    
    
#     sd1 <- sd(pred.1[ , class.lvl ])
#     sd2 <- sd(pred.2[ , class.lvl ])
#     
#     n1 <- length(pred.1[ , class.lvl ])
#     n2 <- length(pred.2[ , class.lvl ])
#     
#     SE <- sqrt( sd1*sd1/n1 + sd2*sd2/n2 )
#     df <- floor( ( sd1^2/n1 + sd2^2/n2 )^2 / 
#                    (((sd1^2 / n1)^2)/(n1 - 1) + ((sd2^2 / n2)^2)/(n2 - 1))
#                  )
#     
#     t.stat <- size/SE
#     p.val <- (1-pt( abs(t.stat) , df=df ))
#     
#     result <- data.frame(var=effect.var, 
#                          value= x, 
#                          size=size, 
#                          t.stat=t.stat, 
#                          p.val=p.val
#                          )
    
    
    return(result)
  })
  
  result <- do.call(rbind, result)
  return(result)
})




