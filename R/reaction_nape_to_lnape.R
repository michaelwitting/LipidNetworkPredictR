#'
#'
#' @export
nape_to_lnape <- function(NAPEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_nape + M_h2o <=> M_lnape + M_fatacid",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  NAPEs <- NAPEs[stringr::str_detect(NAPEs, constraints[1], negate = negate[1])]
  
  
  
  LNAPEs <- unlist(lapply(NAPEs,
                        function(x) {
                          paste0("NAPE(",
                                 lipidomicsUtils::isolate_fatty_acyls(x)[1],
                                 "/0:0/",
                                 lipidomicsUtils::isolate_fatty_acyls(x)[3],
                                 ")")
                        }))
  
  FAs <- unlist(lapply(NAPEs,
                       function(x) {
                         paste0("FA(",
                                lipidomicsUtils::isolate_fatty_acyls(x)[2],
                                ")")
                       }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(NAPEs)),
                            ReactionFormula = character(length(NAPEs)),
                            IsReversible = character(length(NAPEs)),
                            GeneAssociation = character(length(NAPEs)),
                            Pathway = character(length(NAPEs)))
  
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
                                                          "M_nape",
                                                          NAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_lnape",
                                                          LNAPEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          FAs)
  
  
  # return results
  return(list(unique(LNAPEs), unique(FAs), reaction_df))
}