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
        "RHEA:18765", ## cer_to_sm
        "coa_to_acdhap", 
        "RHEA:52716", ## coa_to_fatoh
        "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328", ## coa_to_lpa
        "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715", ## dg_to_sn1mg
        "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278", ## dg_to_sn2mg 
        "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275", ## dg_to_pa
        "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942", ## dg_to_pc
        "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946", ## dg_to_pe
        "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871", ## dg_to_tg
        "RHEA:36179", ## dgo_to_pco
        "RHEA:36187", ## dgo_to_peo
        "dhcer_to_cer",
        "RHEA:44620", ## dhcer_to_dhsm
        "RHEA:19253", ## dhsm_to_dhcer
        "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424", ## fa_to_coa
        "RHEA:45420", ## lnape_to_gpnae
        "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712", ## lpa_to_pa
        "RHEA:36235", ## lpao_to_pao
        "RHEA:15177", ## sn1lpc_to_fa
        "RHEA:44696", ## sn2lpc_to_fa
        "RHEA:12937", ## sn1lpc_to_pc
        "RHEA:32967", ## sn1lpe_to_fa
        "sn2lpe_to_fa", 
        "RHEA:32995", ## sn1lpe_to_pe
        "lpeo_to_peo", 
        "lpep_to_pep", 
        "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466", ## sn1mg_to_dg 
        "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950", ## sn2mg_to_dg 
        "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022", ## sn1mg_to_fa 
        "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874", ## sn2mg_to_fa 
        "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750", ## sn1mg_to_lpa 
        "sn2mg_to_sn1mg", 
        "nae_to_fa", 
        "nape_to_lnape", 
        "nape_to_nae", 
        "nape_to_pnae",
        "napeo_to_nae", 
        "pa_to_cdpdg", 
        "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432", ## pa_to_dg 
        "pao_to_dgo", 
        "RHEA:10604", ## pc_to_dg
        "RHEA:15801", ## pc_to_sn1lpc
        "RHEA:18689", ## pc_to_sn2lpc
        "RHEA:14445", ## pc_to_pa
        "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091", ## pc_to_ps
        "pco_to_lpco", 
        "pe_to_dg", 
        "pe_to_sn1lpe", 
        "RHEA:44408", ## pe_to_sn2lpe
        "pe_to_nape_sn1", 
        "pe_to_nape_sn2", 
        "pe_to_pa", 
        "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609", ## pe_to_ps
        "peo_to_lpeo", 
        "peo_to_napeo_sn1", 
        "peo_to_napeo_sn2", 
        "peo_to_pep", 
        "pep_to_lpep", 
        "pep_to_napep_sn1", 
        "pep_to_napep_sn2", 
        "pg_to_cl", 
        "pgp_to_pg",
        "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831", ## "ps_to_pe", 
        "sm_to_cer", 
        "sphinga_to_dhcer", 
        "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274") ## tg_to_dg)
    
    if (!reaction %in% reaction_types)
        stop(sprintf("%s is not an implemented reaction type", reaction))
    
}