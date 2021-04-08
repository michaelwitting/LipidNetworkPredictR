#'
#'
#' @export
pgp_to_pg <- function(PGPs, template = NA, constraints = c(""), negate = c(FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_h2o + M_pgp <=> M_pi + M_pg",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 1) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PGPs <- PGPs[stringr::str_detect(PGPs, constraints, negate = negate)]
  
  # make new list for LPAs
  PGs <- stringr::str_replace_all(PGPs, "PGP", "PG")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(length(PGs)),
                            ReactionFormula = character(length(PGs)),
                            IsReversible = character(length(PGs)),
                            GeneAssociation = character(length(PGs)),
                            Pathway = character(length(PGs)))
  
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
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pi",
                                                          "Pi")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pgp",
                                                          PGPs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_PG",
                                                          PGs)
  
  # return results
  return(list(unique(PGs), reaction_df))
  
}