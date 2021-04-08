#'
#'
#' @export
lnape_to_gpnae <- function(LNAPEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_lnape + M_h2o <=> M_gpnae + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  LNAPEs <- LNAPEs[stringr::str_detect(LNAPEs, constraints[1], negate = negate[1])]
  
  GPNAEs <- unlist(lapply(LNAPEs,
                          function(x) {
                            paste0("GPNAE(",
                                   lipidomicsUtils::isolate_fatty_acyls(x)[3],
                                   ")")
                          }))
  
  FAs <- unlist(lapply(LNAPEs,
                       function(x) {
                         paste0("FA(",
                                lipidomicsUtils::isolate_fatty_acyls(x)[1],
                                ")")
                       }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(LNAPEs)),
                            ReactionFormula = character(length(LNAPEs)),
                            IsReversible = character(length(LNAPEs)),
                            GeneAssociation = character(length(LNAPEs)),
                            Pathway = character(length(LNAPEs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h2o",
                                                          "H2O")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_lnape",
                                                          LNAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_gpnae",
                                                          LNAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  
  # return results
  return(list(unique(GPNAEs), unique(FAs), reaction_df))
}