#The relevant function in this file is auc(). The other functions are inputs to auc().
#It takes the following arguments

		# value -- a string of the true level of the classification variable we are trying to predict
		# classifications -- a vector of true values in our classification variable
		# probabilities -- a vector of probabilities for our prediction
		
# Example:
	# delta_2x<-auc('Dehisced', bootstrap$Outcome_fin, bootstrap$prob.dehisced.delta3bins2X)




#trapazoid function for calculating the area under a curve
trapezoid = function(x,y){
  len = length(x)
  w <- x[2:len] - x[1:(len-1)]
  h <- (y[2:len] + y[1:(len-1)])/2
  sum(h*w)
}

#Returns Sensitivity and Specificity as a list
sensSpec = function(posValues, probabilities, threshold){
  predValues <- probabilities >= threshold
  sens <- sum(predValues[posValues])/sum(posValues)
  spec <- sum(!predValues[!posValues])/sum(!posValues)
  list(sensitivity=sens, specificity=spec)
}

#Generates ROC curves and AUC
auc = function(value, classifications, probabilities, showPlot=FALSE){
  thresholds <- rev(sort(unique(c(0,probabilities,1))))
  posValues <- classifications==value
  xOneMinusSpec <- vector()
  ySens <- vector()
  i = 1
  for (t in thresholds){
    ss <- sensSpec(posValues, probabilities, t)
    xOneMinusSpec[i] <- 1-ss$specificity
    ySens[i] <- ss$sensitivity
    i = i + 1
  }
  area = trapezoid(xOneMinusSpec, ySens)
  if (showPlot){
    plot(xOneMinusSpec, ySens, type="o", main=paste(value,as.character(signif(area,6)),sep=" - "),
       xlab="1-Specificity", ylab="Sensitivity")
  }
  list(area, data.frame(cbind(xOneMinusSpec,ySens,thresholds)))
}

#source(file="auc.R")