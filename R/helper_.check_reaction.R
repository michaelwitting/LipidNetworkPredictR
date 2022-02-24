#' @name .check_reaction
#' 
#' @title Check if reaction is implemented
#' 
#' @description
#' The helper function checks if the \code{reaction} is in the set of 
#' implemented methods.
#' 
#' The function will throw an error if \code{reaction} is not of length 1 and
#' is not an implemented method.
#' 
#' @details 
#' \code{.check_reaction} is a helper function to test the integrity
#' of \code{reaction}.
#' 
#' @param reaction \code{character(1)}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples
#' .check_reaction("acdhap_to_alkyldhap)
.check_reaction <- function(reaction) {
    
    if (length(reaction) != 1)
        stop("'reaction' has to be of length 1")
    
    ## create vector with the implemented reaction types
    reaction_types <- c("acdhap_to_alkyldhap", "alkyldhap_to_lpao", 
        "c1p_to_cer", "cdpdg_to_pgp", "cdpdg_to_pi", "cer_to_c1p", 
        "cer_to_glccer", "cer_to_sm", "coa_to_acdhap", "coa_to_fatoh",
        "coa_to_lpa", "dg_to_mg", "dg_to_pa", "dg_to_pc", "dg_to_pe",
        "dg_to_tg", "dgo_to_pco", "dgo_to_peo", "dhcer_to_cer", "dhcer_to_dhsm",
        "dhsm_to_dhcer", "fa_to_coa",  "lnapes_to_gpnae", "lpa_to_pa",
        "lpao_to_pao", "sn1lpc_to_fa", "sn21pc_to_fa", "sn1lpc_to_pc", 
        "sn1lpe_to_fa", "sn2lpe_to_fa", "sn1lpe_to_pe",
        "lpeo_to_peo", "lpep_to_pep", "mg_to_dg", "mg_to_fa", "mg_to_lpa",
        "mg_to_mg", "nae_to_fa", "nape_to_lnape", "nape_to_nae", "nape_to_pnae",
        "napeo_to_nae", "pa_to_cdpdg", "pa_to_dg", "pao_to_dgo", "pc_to_dg",
        "pc_to_pa", "pc_to_ps", "pco_to_lpco", "pe_to_dg", "pe_to_sn1lpe",
        "pe_to_nape_sn1", "pe_to_pa", "pe_to_ps", "peo_to_lpeo", "peo_to_napeo_sn1",
        "peo_to_pep", "pep_to_lpep", "pep_to_napep_sn1", "pg_to_cl", "pgp_to_pg",
        "ps_to_pe", "sm_to_cer", "sphinga_to_dhcer", "tg_to_dg")
    
    if (!reaction %in% reaction_types)
        stop(sprintf("%s is not an implemented reaction type", reaction))
    
}