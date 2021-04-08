#'
#'
#' @export
cer_to_glccer <- function(CERs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_udpg + M_c17isocrm <=> M_h + M_udp + M_c17isogluside",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  CERs <- CERs[stringr::str_detect(CERs, constraints, negate = negate)]
  
  # make new list for LCERs
  GLCCERs <- stringr::str_replace_all(CERs, "Cer", "GlcCer")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(GLCCERs)),
                            ReactionFormula = character(length(GLCCERs)),
                            IsReversible = character(length(GLCCERs)),
                            GeneAssociation = character(length(GLCCERs)),
                            Pathway = character(length(GLCCERs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_udpg",
                                                          "UDP-Glucose")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_udp",
                                                          "UDP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isocrm",
                                                          CERs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_c17isogluside",
                                                          GLCCERs)
  
  # return results
  return(list(unique(GLCCERs), reaction_df))
  
}