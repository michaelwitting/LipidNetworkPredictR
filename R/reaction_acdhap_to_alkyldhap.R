#'
#'
#' @export
acdhap_to_alkyldhap <- function(ACDHAPs, FATOHs, template = NA, constraints = c("", ""), negate = c(FALSE, FALSE)) {
  
  # make hard coded template for reaction
  if(is.na(template)) {
    template <- list(reaction_name = "",
                     reaction_formula = "M_adhap + M_alkylR1oh <=> M_h + M_fatacid + M_akdhap",
                     reaction_isReversible = "",
                     reaction_geneAssociation = "",
                     reaction_pathway = "")
  }
  
  # check if constraints have correct lengths
  if(!length(constraints) == 2) {
    stop("Wrong length for constraints")
  }
  
  # constraints variable substrates
  ACDHAPs <- ACDHAPs[stringr::str_detect(ACDHAPs, constraints[1], negate = negate[1])]
  FATOHs <- FATOHs[stringr::str_detect(FATOHs, constraints[2], negate = negate[2])]
  
  # make combinatorics
  ALKYLDHAPs <- data.frame(expand.grid(ACDHAP = ACDHAPs,
                                       FATOH = FATOHs, stringsAsFactors = FALSE),
                           stringsAsFactors = FALSE)
  
  # create new 
  ALKYLDHAPs$FA <- stringr::str_replace_all(ALKYLDHAPs$ACDHAP, "DHAP\\(", "FA(")
  ALKYLDHAPs$ALKYLDHAP <- stringr::str_replace_all(ALKYLDHAPs$FATOH, "FATOH\\(", "DHAP(O-")
  
  # make new data frame with reaction template
  reaction_df <- data.frame(Name = character(nrow(ALKYLDHAPs)),
                            ReactionFormula = character(nrow(ALKYLDHAPs)),
                            IsReversible = character(nrow(ALKYLDHAPs)),
                            GeneAssociation = character(nrow(ALKYLDHAPs)),
                            Pathway = character(nrow(ALKYLDHAPs)))
  
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
  
  
  # adjust variable substrates and products
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_adhap",
                                                          ALKYLDHAPs$ACDHAP)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_alkylR1oh",
                                                          ALKYLDHAPs$FATOH)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_fatacid",
                                                          ALKYLDHAPs$FA)
  
  reaction_df$ReactionFormula <- stringr::str_replace_all(reaction_df$ReactionFormula,
                                                          "M_akdhap",
                                                          ALKYLDHAPs$ALKYLDHAP)
  
  # return results
  return(list(unique(ALKYLDHAPs$ALKYLDHAP), reaction_df))
}