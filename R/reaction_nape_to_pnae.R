#'
#'
#' @export
nape_to_pnae <- function(NAPEs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_nape + M_h2o <=> M_pnae + M_12dag",
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
  
  
  
  PNAEs <- unlist(lapply(NAPEs,
                        function(x) {
                          paste0("PNAE(",
                                 lipidomicsUtils::isolate_fatty_acyls(x)[3],
                                 ")")
                        }))
  
  DGs <- unlist(lapply(NAPEs,
                       function(x) {
                         paste0("DG(",
                                lipidomicsUtils::isolate_fatty_acyls(x)[1],
                                "/",
                                lipidomicsUtils::isolate_fatty_acyls(x)[2],
                                "/0:0)")
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
                                                          "M_pnae",
                                                          PNAEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_12dag",
                                                          DGs)
  
  
  # return results
  return(list(unique(NAEs), unique(DGs), reaction_df))
}