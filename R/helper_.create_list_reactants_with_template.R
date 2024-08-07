#' @name .create_list_products_df_with_template
#' 
#' @title Create the object to return in \code{create_reaction}
#' 
#' @description 
#' Helper function for \code{create_reaction}.
#' 
#' The function \code{.create_list_reactants_with_template} will return a
#' list of length 2 containing the reactants (in the first entry) and the 
#' updated \code{template} object (in the second entry).
#' 
#' @details
#' Depending on the \code{reaction}, the first list entry will contain one or 
#' multiple entries.
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @param df_reaction \code{data.frame}
#' @param template \code{data.frame}
#' 
#' @return list of length 2
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' reaction <- "RHEA:15421"
#' 
#' ## create template
#' template <- LipidNetworkPredictR:::.create_template(template = NA, 
#'     reaction = reaction)
#' 
#' ## create data.frame of substrates
#' df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#'     
#' ## add products to data.frame
#' df_reaction <- LipidNetworkPredictR:::.add_products(
#'     substrates = df_substrates, 
#'     reaction = reaction)
#'     
#' ## make new data.frame with reaction template 
#' LipidNetworkPredictR:::.create_list_reactants_with_template(
#'     df_reaction = df_reaction,
#'     template = template)
.create_list_reactants_with_template <- function(df_reaction, template) {
    
    ## make the entries in df_reaction unique
    l <- lapply(df_reaction, unique)
    
    ## return the list together with the template object as a list
    list(l, template)
    
}