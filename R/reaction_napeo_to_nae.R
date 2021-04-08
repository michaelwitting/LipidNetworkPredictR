#'
#'
#' @export
napeo_to_nae <- function(NAPEOs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_akac2nape + M_h2o <=> M_nae + M_akac2gp",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  NAPEOs <- NAPEOs[stringr::str_detect(NAPEOs, constraints[1], negate = negate[1])]
  
  
  
  NAEs <- unlist(lapply(NAPEOs,
                        function(x) {
                          paste0("NAE(",
                                 lipidomicsUtils::isolate_fatty_acyls(x)[3],
                                 ")")
                        }))
  
  PAOs <- unlist(lapply(NAPEOs,
                       function(x) {
                         paste0("PA(",
                                lipidomicsUtils::isolate_fatty_acyls(x)[1],
                                "/",
                                lipidomicsUtils::isolate_fatty_acyls(x)[2],
                                ")")
                       }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(NAPEOs)),
                            ReactionFormula = character(length(NAPEOs)),
                            IsReversible = character(length(NAPEOs)),
                            GeneAssociation = character(length(NAPEOs)),
                            Pathway = character(length(NAPEOs)))
  
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
                                                          "M_akac2nape",
                                                          NAPEOs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_nae",
                                                          NAEs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akac2gp",
                                                          PAOs)
  
  
  # return results
  return(list(unique(NAEs), unique(PAOs), reaction_df))
}