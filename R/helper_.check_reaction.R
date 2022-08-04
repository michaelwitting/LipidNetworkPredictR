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
#' LipidNetworkPredictR:::.check_reaction(reaction = "RHEA:15421")
.check_reaction <- function(reaction = "RHEA:15421") {
    
    if (length(reaction) != 1)
        stop("'reaction' has to be of length 1")
    
    ## create vector with the implemented reaction types
    reaction_types <- c(
        "RHEA:36171", ## acdhap_to_alkyldhap
        "RHEA:36175", ## alkyldhap_to_lpao
        "c1p_to_cer", 
        "RHEA:12593", ## cdpdg_to_pgp
        "RHEA:11580", ## cdpdg_to_pi
        "RHEA:17929", ## cer_to_c1p
        "RHEA:12088", ## cer_to_glccer
        "cer_to_sm", ## 
        "coa_to_acdhap", 
        "RHEA:52716", ## coa_to_fatoh
        "RHEA:15325", ## coa_to_lpa
        "dg_to_sn1mg", 
        "dg_to_sn2mg", 
        "RHEA:10272", ## dg_to_pa
        "RHEA:32939", ## dg_to_pc
        "RHEA:32943", ## dg_to_pe
        "RHEA:10868", ## dg_to_tg
        "RHEA:36179", ## dgo_to_pco
        "RHEA:36187", ## dgo_to_peo
        "dhcer_to_cer",
        "dhcer_to_dhsm", ## 
        "RHEA:19253", ## dhsm_to_dhcer
        "RHEA:15421", ## fa_to_coa
        "RHEA:45420", ## lnape_to_gpnae
        "RHEA:19709", ## lpa_to_pa
        "RHEA:36235", ## lpao_to_pao
        "RHEA:15177", ## sn1lpc_to_fa
        "RHEA:44696", ## sn2lpc_to_fa
        "RHEA:12937", ## sn1lpc_to_pc
        "RHEA:32967", ## sn1lpe_to_fa
        "sn2lpe_to_fa", 
        "RHEA:32995", ## sn1lpe_to_pe
        "lpeo_to_peo", 
        "lpep_to_pep", 
        "sn1mg_to_dg", 
        "sn2mg_to_dg", 
        "sn1mg_to_fa", 
        "sn2mg_to_fa", 
        "sn1mg_to_lpa", 
        "sn2mg_to_sn1mg", 
        "nae_to_fa", 
        "nape_to_lnape", 
        "nape_to_nae", 
        "nape_to_pnae",
        "napeo_to_nae", 
        "pa_to_cdpdg", 
        "pa_to_dg", 
        "pao_to_dgo", 
        "RHEA:10604", ## pc_to_dg
        "RHEA:15801", ## pc_to_sn1lpc
        "RHEA:18689", ## pc_to_sn2lpc
        "RHEA:14445", ## pc_to_pa
        "RHEA:45088", ## pc_to_ps
        "pco_to_lpco", 
        "pe_to_dg", 
        "pe_to_sn1lpe", 
        "RHEA:44408", ## pe_to_sn2lpe
        "pe_to_nape_sn1", 
        "pe_to_nape_sn2", 
        "pe_to_pa", 
        "RHEA:27606", ## pe_to_ps
        "peo_to_lpeo", 
        "peo_to_napeo_sn1", 
        "peo_to_napeo_sn2", 
        "peo_to_pep", 
        "pep_to_lpep", 
        "pep_to_napep_sn1", 
        "pep_to_napep_sn2", 
        "pg_to_cl", 
        "pgp_to_pg",
        "ps_to_pe", 
        "sm_to_cer", 
        "sphinga_to_dhcer", 
        "tg_to_dg")
    
    if (!reaction %in% reaction_types)
        stop(sprintf("%s is not an implemented reaction type", reaction))
    
}