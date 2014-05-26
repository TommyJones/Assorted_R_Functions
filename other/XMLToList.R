XMLtoList <- function(infilepath){
    ###################################################################################
    # Takes a path to an XML file and turns it into a recursive list of lists
    # in the same format as the XML tree.
    ###################################################################################

    require(XML)
    
    doc <- xmlTreeParse(infilepath, getDTD=FALSE) # read xml into R
    
    r <- xmlRoot(doc) # Gets to the top node in the xml
    
    l <- xmlApply(r, xmlChildren)
        
    parsed <- lapply(l, function(x){
        subparsed <- lapply(x, function(node){
            
            size <- xmlSize(node)
            
            if( size == 0 ){ # If this is an empty node, return NA
                result <- NA        
            }else if(size == 1){ # If this is a "leaf" (no children), return its value
                result <- xmlValue(node)        
            }else{ # Get its children and return their values
                result <- xmlChildren(node)
                result <- lapply(result, xmlValue)
            }
            
            return(result)
        })
        
        return(subparsed)
    })
	
	return(parsed)
}    
