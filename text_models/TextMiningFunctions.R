# Suite of functions to facilitate basic text mining at STPI
# This script focuses on basic text curation and some text clustering

# last update 02/11/2014 - twjones


# Dependent packages
required <- c("tm", "slam", "Matrix", "proxy", "cluster", "compiler", "snowfall")

not.avail <- ! required %in% installed.packages()

if( sum(not.avail) > 0){ install.packages( required ) }

required <- lapply( required , function(x) require(x, character.only=TRUE) )

rm(required, not.avail)

# converts all pdf files in a directory to text files
Pdf2Text <- function(datapath, print=FALSE){
    
    convertToText <- function(name) {
        if(print){ print(name) }
        shell(paste("pdftotext \"", name, "\"", sep=""))
    }
    
    dd <- list.files(datapath, pattern="\\.pdf$", recursive=TRUE, ignore.case=TRUE)
    lapply(dd, function(x) convertToText( paste(datapath, x, sep="") ) )
}



# remove leading/trailing/multiple spaces
FixSpaces <- function(char.vec){
    char.vec <- gsub( "^ +", "", char.vec ) #leading
    char.vec <- gsub( " +$", "", char.vec ) #trailing
    char.vec <- gsub( " +", " ", char.vec ) #multiple to single
    return(char.vec)
}

# makes some adjustments to pluralization for a vector whose entries are single words (not whole documents).
# WARNING: This does make mistakes for irregular words. You should check its results manually.
CorrectS <- function(term.vec){
    s.adjust <- gsub("sses$", "ss", term.vec) 
    keep.list <- s.adjust[ grepl("sis$|ss$|us$", s.adjust) | nchar(term.vec) <= 3 | grepl( "_[a-zA-Z][a-zA-Z][a-zA-Z]$", s.adjust) ]
    
    s.adjust2 <- gsub("ies$", "y", s.adjust)
    s.adjust2 <- gsub("s$", "", s.adjust2)
    
    out.list <- s.adjust2
    out.list[ s.adjust %in% keep.list ] <- s.adjust[ s.adjust %in% keep.list ]
    
    result <- data.frame(original=term.vec, adjusted=out.list, changed=term.vec!=out.list, stringsAsFactors=FALSE)
    return(result)
}


# Make a document term matrix and get term frequencies
TermFreqDTM <- function(corpus){
    # input is a charactar vector where each entry is a document with tokens separated by a single space
    # output is a list with two entries
    # term.freq is a data frame with words, term frequencies, and document frequencies
    # dtm is a document term matrix in simple triplet form as given by the tm package
    
    require(tm)
    require(slam)
    
    tmcorp <- Corpus(VectorSource(corpus))
    tmdtm <- DocumentTermMatrix(tmcorp)
    
    term.doc.freq <- col_sums( tmdtm > 0 )
    term.freq <- col_sums( tmdtm )
    
    term.freq <- data.frame(term=names(term.freq), term.freq=term.freq, doc.freq=term.doc.freq, stringsAsFactors=FALSE)
    
    result <- list(term.freq=term.freq, dtm=tmdtm)
    
    return(result)
}


MakeSparseDTM <- function(dtm){
    # converts a simple triplet DocumentTermMatrix object from the tm package to a sparse matrix from the Matrix package
    # dtm is a simple triplet matrix
    dtm.sparse <- sparseMatrix(i=dtm$i, j=dtm$j, x=dtm$v, 
                               dims=c(dtm$nrow, dtm$ncol))
    
    rownames(dtm.sparse) <- Docs(dtm)
    colnames(dtm.sparse) <- Terms(dtm)
    
    return(dtm.sparse)
}

GetNN<-cmpfun(function(sim.mat,N, edgelist=TRUE){
    # takes in a symmetric matrix of similarities and returns a list of matrices were all entries in the similarity 
    # matrix are set to 0 except the top N similarites.
    # if edgelist is TRUE (default) a data frame is returned in the form of an edgelist
    # otherwise, a sparse adjacency matrix is returned
    
    diag(sim.mat)<-0
    
    tmp<-apply(sim.mat,2,function(x){
        top.N <-sort(x,decreasing=TRUE)[1:N]
        
        x[!names(x) %in% names(top.N)]<-0
        
        return(x)
    })
    
    #tmp is not symmetric; make it that way. Think of tmp as a directed adjacency matrix. We want it undirected.
    result <- pmax(tmp,t(tmp)) #takes the max of two values. If both 0, then 0. if both non-zero, then they're equal and non-zero. If one is zero and one is non-zero, return the non-zero value.
    
    if(! edgelist ){ 
        return(result) 
    }else{
        require( igraph )
        result.graph <- graph.adjacency(result, mode="undirected", weighted=TRUE)
        result <- as.data.frame( get.edgelist(result.graph, names=TRUE))
        names(result) <- c("Source", "Target")
        result$weight <- E(result.graph)$weight
        return(result)
    }
})

# association score of terms in documents
# terms can be a single character or a vector of characters
TermSearch <- function( terms, dtm.sparse ){
    
    terms <- terms[ terms %in% colnames(dtm.sparse) ]
        
    pmat <- dtm.sparse / rowSums(dtm.sparse)
    
    if( length(terms) > 1 ){ 
        pterm <- sum(rowSums(dtm.sparse[ , terms ])) / sum(colSums(dtm.sparse))
        result <- rowSums(pmat[ , terms ]) - pterm 
    }else{
        pterm <- sum(dtm.sparse[ , terms ]) / sum(colSums(dtm.sparse))
        result <- as.numeric(pmat[ , terms ] - pterm)
    }
    
    names(result) <- rownames(dtm.sparse)
    
    return(result[ order(result, decreasing=TRUE ) ])
}


# returns cosine similarity of two vectors, x and y
# can be used alone or accompanying the 'dist' function from the proxy package, loaded with this file
CosineSim <- cmpfun(function(x,y){
    result <- sum(x*y) / sqrt( sum(x*x) * sum(y*y) )
})




