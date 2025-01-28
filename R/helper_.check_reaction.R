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
        "RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275", ## dg_to_pa
        "RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607", ## pc_to_dg
        "RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871", ## dg_to_tg
        "RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583", ## cdpdg_to_pi
        "RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091", ## cer_to_glccer
        "RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596", ## cdpdg_to_pgp
        "RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940", ## sn1lpc_to_pc
        "RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448", ## pc_to_pa
        "RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180", ## sn1lpc_to_fa
        "RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328", ## coa_to_lpa
        "RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424", ## fa_to_coa
        "RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804", ## pc_to_sn1lpc
        "RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232", ## pa_to_cdpdg
        "RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248", ## lpep_to_pep 
        "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744", ## sn2mg_to_dg
        "RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908", ## lpep_to_fal
        "RHEA:17505", "RHEA:17506", "RHEA:17507", "RHEA:17508", ## nae_to_fa
        "RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660", ## coa_to_acyldhap
        "RHEA:17729", "RHEA:17730", "RHEA:17731", "RHEA:17732", ## coa_to_ce
        "RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932", ## cer_to_cerp
        "RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004", ## pi_to_sn1lpi
        "RHEA:18689", "RHEA:18690", "RHEA:18691", "RHEA:18692", ## pc_to_sn2lpc
        "RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768", ## cer_to_sm
        "RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712", ## lpa_to_pa
        "RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831", ## ps_to_pe 
        "RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959", ## peo_to_pep, 
        "RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995", ## lpco_to_pco
        "RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432", ## pa_to_dg 
        "RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609", ## pe_to_ps
        "RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874", ## sn2mg_to_fa 
        "RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934", ## pg_to_cl
        "RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938", ## cl_to_lcl 
        "RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942", ## dg_to_pc
        "RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946", ## dg_to_pe
        "RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950", ## sn2mg_to_dg
        "RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970", ## sn1lpe_to_fa
        "RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998", ## sn1lpe_to_pe
        "RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162", ## nape_to_nae
        "RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198", ## sn1lpi_to_pi
        "RHEA:33203", "RHEA:33204", "RHEA:33205", "RHEA:33206", ## lpg_to_pg
        "RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274", ## tg_to_dg
        "RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278", ## dg_to_sn2mg 
        "RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746", ## cerp_to_cer 
        "RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750", ## sn1mg_to_lpa 
        "RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754", ## pgp_to_pg,
        "RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022", ## sn1mg_to_fa 
        "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666", ## dg_to_sn1mg
        "RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842", ## lcl_to_cl
        "RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086", ## lpco_to_mgo
        "RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174", ## acyldhap_to_alkyldhap
        "RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178", ## alkyldhap_to_lpao
        "RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182", ## dgo_to_pco
        "RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190", ## dgo_to_peo
        "RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198", ## pep_to_lpep,
        "RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202", ## lpep_to_mgp
        "RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206", ## lpep_to_lpap
        "RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234", ## pco_to_lpco 
        "RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238", ## lpao_to_pao
        "RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242", ## pao_to_dgo 
        "RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466", ## sn1mg_to_dg
        "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886", ## fa_to_coa
        "RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930", ## lpco_to_lpao
        "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946", ## sn1mg_to_dg
        "RHEA:39995", "RHEA:39996", "RHEA:39997", "RHEA:39998", ## nae_to_fa
        "RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487", ## pi_to_dg
        "RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411", ## pe_to_sn2lpe
        "RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607", ## pe_to_sn1lpe 
        "RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623", ## dhcer_to_dhsm
        "RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699", ## sn2lpc_to_fa
        "RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715", ## dg_to_sn1mg
        "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867", ## tg_to_dg
        "RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091", ## pc_to_ps
        "RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191", ## pe_to_nape_sn1
        "RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195", ## pe_to_nape_sn2
        "RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303", ## dhsm_to_dhcer
        "RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423", ## lnape_to_gpnae
        "RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463", ## nape_to_lnape 
        "RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547", ## dhcer_to_cer
        "RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647", ## sm_to_cer 
        "RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719", ## coa_to_fao
        "RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427", ## sphinga_to_dhcer 
        "RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599", ## pep_to_napep_sn1
        "RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954", ## pe_to_dg
        "sn2lpe_to_fa", 
        "lpeo_to_peo", 
        "sn2mg_to_sn1mg", 
        "nape_to_pnae",
        "napeo_to_nae", 
        "pe_to_pa", 
        "peo_to_lpeo", 
        "peo_to_napeo_sn1", 
        "peo_to_napeo_sn2", 
        "pep_to_napep_sn2")
    
    if (!reaction %in% reaction_types)
        stop(sprintf("%s is not an implemented reaction type", reaction))
    
}