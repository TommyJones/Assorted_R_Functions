##############################################
# functions for cross validation
##############################################

CvData <- function(data, k, seed=NULL){
    
    # set seed if desired / not already set
    if( ! is.null(seed) ) set.seed(seed)
    
    # set up basics
    rows <- 1:nrow(data)
    
    size <- floor(nrow(data)/k)
    
    remainder <- nrow(data) - size * k
    
    # get the indexes of each fold into a list
    folds <- vector(mode="list", length=k)
    
    for(j in 1:k){
        folds[[ j ]] <- sample(rows, size=size, replace=FALSE)
        rows <- rows[ ! rows %in% folds[[ j ]] ]
    }
    
    if( remainder > 0){
        for(j in 1:length(rows)){
            folds[[ j ]][ length(folds[[ j ]]) + 1 ] <- sample(rows, size=1)
            rows <- rows[ ! rows %in% folds[[ j ]] ]
        }
    }
    
    # combinations of each fold to use for training sets / test sets
    combos <- combn(x=1:k, m=(k-1))
    
    # pull indices into format that's easier to use
    cv.folds <- vector(mode="list", length=k)
    
    for( j in 1:k ){
        cv.folds[[ j ]]$training <- unlist( folds[ combos[ , j ] ] )
        cv.folds[[ j ]]$test <- folds[[ (1:k)[ ! (1:k) %in% combos[ , j ] ]  ]]
    }
    
    result <- list(fold.indices=cv.folds, data=data)
    return(result)
}








