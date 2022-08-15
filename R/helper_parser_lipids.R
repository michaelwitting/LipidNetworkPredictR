#' @title Get all fatty acyls, alkyls and alkenyls
#' 
#' @description
#' This functions isolates all fatty acyls, alkyls and alkenyls from a given 
#' lipid shorthand notation and returns them as vector. 
#' 
#' @details
#' Supported modifications are currently hydroxy groups (OH), hydroperoxy 
#' groups (OOH), keto groups (O), and amino groups (NH2).
#' 
#' @param lipids list of vector of shorthand notation of a acyl (as string), 
#' e.g. \code{c("PC(18:0/20:4(7E,9E,11Z,14Z))", "PC(16:0/18:1(9Z))")}
#' 
#' @examples
#' lipids <- c("PC(18:0/20:4(7E,9E,11Z,14Z))", "PC(16:0/18:1(9Z))")
#' isolate_radyls(lipids)
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @importFrom stringi stri_extract_all_regex stri_detect_regex
#' @export
isolate_radyls <- function(lipids) {
    
    ## iterate over lipds and isolate radyls from single lipids, 
    ## use of lapply for list and vectors
    lapply(lipids, FUN = .isolate_radyls_helper)
}

## helper function to isolate radyls from a single lipid
.isolate_radyls_helper <- function(lipid) {
    
    ## get all possible building blocks
    radyls <- stringi::stri_extract_all_regex(str = lipid, 
        pattern = "(m|d|t|O-|P-)*\\d+:\\d+(\\((\\d*(E|Z|Me|OH|OOH|O|NH2|delta)(\\[(S|R)\\])*,*)*\\))*")[[1]]
    
    ## remove sphingoid bases
    filter <- stringi::stri_detect_regex(str = radyls, pattern = "(m|d|t)", 
        negate = TRUE)
    radyls <- radyls[filter]
    
    ## return result
    return(radyls)
}