#'
#'
#' @export
pg_to_cl <- function(PGs, CDPDGs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_cdpdag + M_pg <=> M_h + M_cmp + M_clpn",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  PGs <- PGs[stringr::str_detect(PGs, constraints[1], negate = negate[1])]
  CDPDGs <- CDPDGs[stringr::str_detect(CDPDGs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  CLs <- data.frame(expand.grid(PGs = PGs,
                                CDPDGs = CDPDGs, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  
  # isolate core from PGs
  CLs$PGs2 <- stringr::str_replace_all(CLs$PGs, "^PG\\(", "")
  CLs$PGs2 <- stringr::str_replace_all(CLs$PGs2, "\\)$", "")
  
  #isolate core from CDPDGs
  CLs$CDPDGs2 <- stringr::str_replace_all(CLs$CDPDGs, "^CDPDG\\(", "")
  CLs$CDPDGs2 <- stringr::str_replace_all(CLs$CDPDGs2, "\\)$", "")
  
  #create CL
  CLs$CL <- paste0("CL(1'-[", CLs$PGs2,"],3'-[", CLs$CDPDGs2, "])")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(CLs)),
                            ReactionFormula = character(nrow(CLs)),
                            IsReversible = character(nrow(CLs)),
                            GeneAssociation = character(nrow(CLs)),
                            Pathway = character(nrow(CLs)))
  
  # fill with data
  reaction_df$Name <- template$reaction_name
  reaction_df$ReactionFormula <- template$reaction_formula
  reaction_df$IsReversible <- template$reaction_isReversible
  reaction_df$GeneAssociation <- template$reaction_geneAssociation
  reaction_df$Pathway <- template$reaction_pathway
  
  # adjust fixed substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_h",
                                                          "H+")
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cmp",
                                                          "CMP")
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_cdpdag",
                                                          CLs$CDPDGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_pg",
                                                          CLs$PGs)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_clpn",
                                                          CLs$CL)
  
  # return results
  return(list(unique(CLs$CL), reaction_df))
}