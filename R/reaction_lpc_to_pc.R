#'
#'
#' @export
sn1lpc_to_pc <- function(sn1LPCs, CoAs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_ag3pc + M_fataccoa <=> M_pchol + M_coa",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  sn1LPCs <- sn1LPCs[stringr::str_detect(sn1LPCs, constraints[1], negate = negate[1])]
  CoAs <- CoAs[stringr::str_detect(CoAs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  PCs <- data.frame(expand.grid(sn1LPCs = sn1LPCs,
                                CoAs = CoAs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  
  PCs$PCs <- stringr::str_replace_all(PCs$sn1LPCs,
                                      "/0:0",
                                      unlist(lapply(PCs$CoAs,
                                                    function(x) {
                                                      paste0("/",lipidomicsUtils::isolate_fatty_acyls(x))
                                                    }
                                      )))
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(PCs)),
                            ReactionFormula = character(nrow(PCs)),
                            IsReversible = character(nrow(PCs)),
                            GeneAssociation = character(nrow(PCs)),
                            Pathway = character(nrow(PCs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_coa",
                                                          "CoA")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fataccoa",
                                                          PCs$CoAs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_ag3pc",
                                                          PCs$sn1LPCs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pchol",
                                                          PCs$PCs)
  
  # return results
  return(list(unique(PCs$PCs), reaction_df))
}
