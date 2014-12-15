require(tm)
require(slam)
require(Matrix)
require(snowfall)

options(mc.cores=1)

Files2Vec <- function(directory){
	file.list <- grep("\\.txt", dir(directory))
	
	vec <- sapply(file.list, function(x){
		result <- scan(paste(directory, x, sep=""), what="character", sep="\n")
		result <- gsub("\\s", " ", result)
		return(result)
	})
	
	names(vec) <- gsub("[^a-zA-Z]", "_", file.list)
	
	return(vec)
}

CorrectS <- function(term.vec){
	# makes some adjustments to pluralization
	# WARNING: This does make mistakes for irregular words. You should check its results manually.
    s.adjust <- gsub("sses$", "ss", term.vec) 
    keep.list <- s.adjust[ grepl("sis$|ss$|us$|species$", s.adjust) | nchar(term.vec) <= 3 | grepl( "_[a-zA-Z][a-zA-Z][a-zA-Z]$", s.adjust) ]
    
    s.adjust2 <- gsub("ies$", "y", s.adjust)
    s.adjust2 <- gsub("s$", "", s.adjust2)
    
    out.list <- s.adjust2
    out.list[ s.adjust %in% keep.list ] <- s.adjust[ s.adjust %in% keep.list ]
    
    result <- data.frame(original=term.vec, adjusted=out.list, changed=term.vec!=out.list, stringsAsFactors=FALSE)
    return(result)
}

MakeSparseDTM <- function(dtm){
  # dtm is a simple triplet matrix
  dtm.sparse <- sparseMatrix(i=dtm$i, j=dtm$j, x=dtm$v, 
                             dims=c(dtm$nrow, dtm$ncol))
  
  rownames(dtm.sparse) <- Docs(dtm)
  colnames(dtm.sparse) <- Terms(dtm)
  
  return(dtm.sparse)
}

NgramTokenizer <- function(min, max) {
  require(RWeka)
  # Function creates a function to create ngrams from a document term matrix
  # For bigrams min=max=2. For bigrams and trigrams min=2, max=3
  # Example: 
  # Bigrams <- NgramTokenizer(2, 2)
  # myDTM <- DocumentTermMatrix(myCorp, control = list(tokenize = Bigrams))
  
  
  result <- function(x) {NGramTokenizer(x, Weka_control(min = min, max = max))}
  
  return(result)
}

Vec2Dtm <- function(vec, max.n.gram=1, remove.stopwords=TRUE, lower=TRUE, remove.punctuation=TRUE, remove.numbers=TRUE, stem.document=FALSE, custom.stopwords=NULL){
	# for now, it is strongly advised to accept the defaults for lower, remove.punctuation, and remove.numbers
	# Other functions are built assuming that the column headings of a dtm contain only letters and underscores "_"
	
	stopwords <- unique(c(stopwords("english"), stopwords("SMART")))
    
    if( ! is.null(custom.stopwords) ) stopwords <- custom.stopwords
	
	if( lower ) vec <- tolower(vec)
	
	if( remove.punctuation ){ 
		vec <- gsub("[^a-zA-Z0-9]", " ", vec)
		stopwords <- gsub("[^a-zA-Z0-9]", " ", stopwords)
		stopwords <- unique(unlist(strsplit(stopwords, split="\\s+")))
	}
	
	if( remove.numbers ){ 
		vec <- gsub("[0-9]", " ", vec)
	}

	vec <- gsub("\\s+", " ", vec) # remove extra spaces
	
	corp <- Corpus(VectorSource(vec))
    
    
    if( remove.stopwords ){
        corp <- tm_map(x=corp, removeWords, stopwords)
        corp <- tm_map(x = corp, stripWhitespace)
    }
    
    if(stem.document){
        corp <- tm_map(x=corp, stemDocument)
    }
	
    if( max.n.gram > 1 ){
        dtm <- DocumentTermMatrix(corp, control=list(tokenize=NgramTokenizer(min=1, max=max.n.gram)))
    }else{
        dtm <- DocumentTermMatrix(corp)
    }
	
	dtm <- MakeSparseDTM(dtm=dtm)
	
	colnames(dtm) <- gsub(" ", "_", colnames(dtm))
	
	return(dtm)
}

DepluralizeDtm <- function(dtm){

	# run depluralization on column names of dtm
	adjust.terms <- CorrectS(colnames(dtm))
	
	# begin procedure to merge depluralized columns with their singular form

    colnames(dtm) <- adjust.terms$adjusted

	# partition dtm.tmp on columns with duplicates (terms that had plurals) and columns without
	unchanged <- dtm[ , ! colnames(dtm) %in% colnames(dtm)[ duplicated(colnames(dtm)) ] ]
	changed <- dtm[ , colnames(dtm) %in% colnames(dtm)[ duplicated(colnames(dtm)) ] ]

	# get indices of columns to be merged
	sfInit(parallel=TRUE, cpus=4)
	sfExport("changed")
    sfLibrary(Matrix)
	
	term.indices <- sfLapply( unique(colnames(changed)), function(x) grep(x, colnames(changed), fixed=TRUE))
	
    sfStop()
    
    gc()
    
	# merge columns by summation
	sfInit(parallel=TRUE, cpus=4)
	sfExport("changed")
	sfLibrary(Matrix)
	
	temp <- sfLapply( term.indices, function(x){
		result <- Matrix(rowSums(changed[ , x ]), sparse=TRUE)
	})
    
    sfStop()

    gc()

	# put back together into a sparse matrix
	batches <- seq(1, length(temp), by=100)
    
    sfInit(parallel=TRUE, cpus=4)
    sfExport("temp")
	sfLibrary(Matrix)

	temp2 <- sfLapply(batches, function(x){
		do.call(cBind, temp[ x:min( x + 99, length(temp) ) ])
	})
    
    sfStop()

    if(length(temp2) > 100 ){
        
        batches <- seq(1, length(temp2), by=100)
        
        sfInit(parallel=TRUE, cpus=4)
        sfExport("temp2")
        sfLibrary(Matrix)
        
        temp <- sfLapply(batches, function(x){
            do.call(cBind, temp2[ x:min( x + 99, length(temp2) ) ])
        })
        
        sfStop()
        
        gc()
        
        temp <- do.call(cBind, temp)
        
    }else{
        temp <- do.call(cBind, temp2)
    }

    colnames(temp) <- unique(colnames(changed))

	dtm <- cBind(unchanged, temp)

	return(dtm)
}

TermDocFreq <- function(dtm){
	freq.mat <- data.frame(term=colnames(dtm), term.freq=colSums(dtm), doc.freq=colSums(dtm > 0), stringsAsFactors=FALSE)
	freq.mat$idf <- log(nrow(dtm) / freq.mat$doc.freq)
	return(freq.mat)
}

CountTerms <- function(freq.mat, max.docs, min.docs){
	keep.terms <- freq.mat$term[ freq.mat$doc.freq > min.docs & freq.mat$doc.freq < max.docs ]
	keep.index <- freq.mat$term %in% keep.terms
	result <- list(keep.terms=keep.terms, 
					summary.tf=summary(freq.mat$term.freq[ keep.index ]),
					summary.df=summary(freq.mat$doc.freq[ keep.index ]),
					summary.idf=summary(freq.mat$idf[ keep.index ]),
					num.terms=length(keep.terms), 
					num.uni=sum(! grepl("_", keep.terms)),
					num.bi=sum( grepl("^[a-z]+_[a-z]+$", keep.terms)),
					num.tri=sum( grepl("^[a-z]+_[a-z]+_[a-z]+$", keep.terms)),
					num.more=sum( grepl("^[a-z]+_[a-z]+_[a-z]+_[a-z]", keep.terms) )
					)
	return(result)
}

SubsetDtm <- function(dtm, keep.terms){
	dtm <- dtm[ , keep.terms ]
	return(dtm)
}
