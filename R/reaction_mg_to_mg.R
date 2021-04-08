#'
#'
#' @export
sn2mg_to_sn1mg <- function(sn2MGs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_mag <=> M_1magol",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn2MGs <- sn2MGs[stringr::str_detect(sn2MGs, constraints[1], negate = negate[1])]
  
  # "sn2 loss"
  sn1MGs <- unlist(lapply(sn2MGs, function(sn2mg) {
    fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(sn2mg)
    fa <- paste0("MG(", fatty_acyls[2], "/0:0/0:0)")
    return(fa)
  }))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(sn2MGs)),
                            ReactionFormula = character(length(sn2MGs)),
                            IsReversible = character(length(sn2MGs)),
                            GeneAssociation = character(length(sn2MGs)),
                            Pathway = character(length(sn2MGs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway

  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_mag",
                                                          sn2MGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_1magol",
                                                          sn1MGs)
  
  # return results
  return(list(unique(sn1MGs),reaction_df))
}
