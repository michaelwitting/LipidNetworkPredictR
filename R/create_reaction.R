#' @name create_reaction
#' 
#' @title Create reaction for a given reaction of lipid metabolism
#' 
#' @description 
#' The function \code{create_reaction} creates from a given template reaction
#' the specific reaction given a set of substrates.
#' 
#' The function returns a list of the products of the reaction (in the first
#' entry of the returned list), together with the \code{template} with the 
#' updated reactions (in the second entry of the returned list). 
#' 
#' @details 
#' If \code{template == NULL}, the function will load a template specific
#' to the \code{reaction} with mininum information on the reaction.
#' 
#' @param substrates \code{list}, with entry names equal to substrates
#' @param template \code{data.frame}
#' @param reaction \code{character(1)}
#' @param constraints \code{character}, with length \code{length(substrates)}
#' @param negate \code{logical}, with length \code{length(substrates)}
#' 
#' @export
#' @return list of length 2
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'    and Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#'     
#' template <- list(reaction_name = "CoA biosynthesis", 
#'     reaction_formula = "M_atp + M_coa + M_fatacid <=> M_ppi + M_amp + M_fataccoa",
#'     reaction_isReversible = "",
#'     reaction_geneAssociation = "",
#'     reaction_pathway = "Acyl-CoA synthetase")
#'  
#' ## run create_reaction (template = NULL)
#' create_reaction(substrates = list(FA = FA), template = NULL,
#'     reaction = "fa_to_coa")
#'     
#' ## run create_reaction (template = template)
#' create_reaction(substrates = list(FA = FA), template = template,
#'     reaction = "fa_to_coa")
create_reaction <- function(substrates, template = NULL, reaction = "fa_to_coa", 
    constraints = character(length(substrates)), 
    negate = logical(length(substrates))) {
    
    ## create template
    template <- .create_template(template = template, reaction = reaction)
    
    ## create data.frame of substrates
    df_substrates <- .create_substrates_combinations(substrates = substrates,
        constraints = constraints, negate = negate)
    .check_colnames_substrates_combinations(df = df_substrates,
        reaction = reaction)
    
    ## add products to data.frame
    df_reaction <- .add_products(substrates = df_substrates, reaction = reaction)
    
    ## make new data.frame with reaction template
    df <- .create_df_with_template(df_reaction = df_reaction, 
        template = template, reaction = reaction)
    
    ## return results
    .create_list_products_df_with_template(df_reaction = df_reaction, 
        template = df, reaction = reaction)
}
