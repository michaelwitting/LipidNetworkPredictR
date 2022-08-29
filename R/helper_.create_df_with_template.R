
#' @name .create_df_with_template
#' 
#' @title Create data frame and substitute substrates/products
#' 
#' @description 
#' Helper function for \code{create_reaction}.
#' 
#' The function \code{.create_df_with_template} takes the \code{template}
#' and replaces the \code{"reactionFormula"} by the substrates/products from
#' \code{df_reaction}.
#' 
#' @details
#' The string replacement depend on the \code{reaction} argument.
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de} and 
#'     Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @param df_reaction \code{data.frame}
#' @param template \code{df_template}
#' @param reaction \code{character(1)}
#' 
#' @return data.frame
#' 
#' @importFrom stringi stri_replace_all_fixed stri_replace_first_fixed
#' @importFrom stringi stri_replace_all_regex
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' reaction <- "RHEA:15421"
#' 
#' ## create template
#' template <- LipidNetworkPredictR:::.create_template(template = NA, 
#'     reaction = reaction)
#' 
#' ## create data.frame of substrates
#' df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#'     
#' ## add products to data.frame
#' df_reaction <- LipidNetworkPredictR:::.add_products(
#'     substrates = df_substrates, 
#'     reaction = reaction)
#'     
#' ## make new data.frame with reaction template
#' LipidNetworkPredictR:::.create_df_with_template(
#'     df_reaction = df_reaction,
#'     template = template, reaction = reaction)
.create_df_with_template <- function(df_reaction, template, 
    reaction = "RHEA:15421") {
 
    .check_reaction(reaction)
    
    ## create a data.frame and will with ""
    .char <- character(nrow(df_reaction))
    
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274",
            "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867")) 
        .char <- c(.char, .char)
    
    df <- data.frame(reaction_name = .char, reaction_formula = .char,
        reaction_isReversible = .char, reaction_geneAssociation = .char, 
        reaction_pathway = .char)
    
    ## fill with data
    df[["reaction_name"]] <- template[["reaction_name"]]
    df[["reaction_formula"]] <- template[["reaction_formula"]]
    df[["reaction_isReversible"]] <- template[["reaction_isReversible"]]
    df[["reaction_geneAssociation"]] <- template[["reaction_geneAssociation"]]
    df[["reaction_pathway"]] <- template[["reaction_pathway"]]
    
    ## get the (reaction) formula
    .formula <- df$reaction_formula
    
    ## depending on the reaction, adjust fixed and variable substrates and 
    ## products
    
    ## acdhap_to_alkyldhap
    if (reaction %in% c("RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylDHAP", replacement = df_reaction[["AcylDHAP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FATOH", replacement = df_reaction[["FATOH"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AlkylDHAP", replacement <- df_reaction[["AlkylDHAP"]])
    }
    
    ## alkyldhap_to_lpao
    if (reaction %in% c("RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AlkylDHAP", replacement = df_reaction[["AlkylDHAP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPA-O", replacement = df_reaction[["LPAO"]])
    }
    
    if (reaction == "c1p_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_C1P", replacement = df_reaction[["C1P"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
    }
    
    ## cdpdg_to_pgp
    if (reaction %in% c("RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CDP-DG", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PGP", replacement = df_reaction[["PGP"]])
    }
    
    ## cdpdg_to_pi
    if (reaction %in% c("RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CDP-DG", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PI", replacement = df_reaction[["PI"]])
    }
    
    if (reaction == "RHEA:17929") { ## cer_to_c1p
        .formula <- stringi::stri_replace_first_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_C1P", replacement = df_reaction[["C1P"]])
    }
    
    if (reaction == "RHEA:12088") { ## cer_to_glccer
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_GLCCER", replacement = df_reaction[["GLCCER"]])
    }
    
    if (reaction == "RHEA:18765") { ## cer_to_sm
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_SM", replacement = df_reaction[["SM"]])
    }
    
    ## cl_to_lcl
    if (reaction %in% c("RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CL", replacement = df_reaction[["CL"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2,4-LCL", replacement = df_reaction[["LCL"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## coa_to_acyldhap
    if (reaction %in% c("RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylDHAP", replacement = df_reaction[["AcylDHAP"]])
    }
    
    ## coa_to_fatoh
    if (reaction %in% c("RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FATOH",  replacement = df_reaction[["FATOH"]])
    }
    
    ## coa_to_lpa
    if (reaction %in% c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPA", replacement = df_reaction[["LPA"]])
    }
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715",
            "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
                pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1-MG", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## dg_to_sn2mg
    if (reaction %in% c("RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-MG", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## dg_to_pa
    if (reaction %in% c("RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275")) {  
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PA", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## dg_to_pc
    if (reaction %in% c("RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
    }
    
    ## dg_to_pe
    if (reaction %in% c("RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
    }
    
    ## dg_to_tg
    if (reaction %in% c("RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_TG", replacement = df_reaction[["TG"]])
    }
    
    ## dgo_to_pco
    if (reaction %in% c("RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182")) {
        .formula <- stringi::stri_replace_first_fixed(str = .formula, 
            pattern = "M_DG-O", replacement = df_reaction[["DGO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC-O", replacement = df_reaction[["PCO"]])
    }

    ## dgo_to_peo
    if (reaction %in% c("RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190")) {
        .formula <- stringi::stri_replace_first_fixed(str = .formula,
            pattern = "M_DG-O", replacement = df_reaction[["DGO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
    }
    
    if (reaction == "dhcer_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHCER", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
    }
    
    if (reaction == "RHEA:44620") { ## dhcer_to_dhsm
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHCER", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHSM", replacement = df_reaction[["DHSM"]])
    }
    
    if (reaction == "RHEA:19253") { ## dhsm_to_dhcer
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHCER", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHSM", replacement = df_reaction[["DHSM"]])
    }
    
    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424",
            "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
    }
    
    ## lcl_to_cl
    if (reaction %in% c("RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2,4-LCL", replacement = df_reaction[["LCL"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CL", replacement = df_reaction[["CL"]])
    }
    
    if (reaction == "RHEA:45420") { ## lnape_to_gpnae
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LNAPE", replacement = df_reaction[["LNAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_gpnae", replacement = df_reaction[["GPNAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## lpa_to_pa
    if (reaction %in% c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPA", replacement = df_reaction[["LPA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PA", replacement = df_reaction[["PA"]])
    }
    
    ## lpao_to_pao
    if (reaction %in% c("RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_LPA-O", replacement = df_reaction[["LPAO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA-O", replacement = df_reaction[["PAO"]])
    }
    
    ## sn1lpc_to_fa
    if (reaction %in% c("RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## sn2lpc_to_fa
    if (reaction %in% c("RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPC", replacement = df_reaction[["sn2LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## sn1lpc_to_pc
    if (reaction %in% c("RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
    }
    
    ## sn1lpe_to_fa
    if (reaction %in% c("RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1-LPE", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "sn2lpe_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPE", replacement = df_reaction[["sn2LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## sn1lpe_to_pe
    if (reaction %in% c("RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPE", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
    }
    
    ## sn1lpi_to_pi
    if (reaction %in% c("RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
                pattern = "M_1-LPI", replacement = df_reaction[["sn1LPI"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
                pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
                pattern = "M_PI", replacement = df_reaction[["PI"]])
    }
    
    if (reaction == "lpeo_to_peo") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ak2lgpe", replacement = df_reaction[["LPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
    }
    
    ## lpep_to_pep
    if (reaction %in% c("RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPE-P", replacement = df_reaction[["LPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-P", replacement = df_reaction[["PEP"]])
    }
    
    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466",
            "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-MG", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950",
            "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-MG", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## sn1mg_to_fa
    if (reaction %in% c("RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-MG", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])   
    }
    ## sn2mg_to_fa
    if (reaction %in% c("RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-MG", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## sn1mg_to_lpa
    if (reaction %in% 
            c("RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-MG", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPA", replacement = df_reaction[["LPA"]])
    }
    
    if (reaction == "sn2mg_to_sn1mg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-MG", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-MG", replacement = df_reaction[["sn1MG"]])
    }
    
    if (reaction == "nae_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAE", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "nape_to_lnape") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPE", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LNAPE", replacement = df_reaction[["LNAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "nape_to_nae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_NAPE", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_NAE", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA", replacement = df_reaction[["PA"]])
    }
    
    if (reaction == "nape_to_pnae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_NAPE", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PNAE", replacement = df_reaction[["PNAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "napeo_to_nae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_NAPEO", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_NAE", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA-O", replacement = df_reaction[["PAO"]])
    }
    
    ## pa_to_cdpdg
    if (reaction %in% c("RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_CDP-DG", replacement = df_reaction[["CDPDG"]])
    }
    
    ## pa_to_dg
    if (reaction %in% c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## pao_to_dgo
    if (reaction %in% c("RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PA-O", replacement = df_reaction[["PAO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DG-O", replacement = df_reaction[["DGO"]])
    }
    
    ## pc_to_dg
    if (reaction %in% c("RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## pc_to_sn1lpc
    if (reaction %in% c("RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804")) { 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## pc_to_sn2lpc
    if (reaction %in% c("RHEA:18689", "RHEA:18690","RHEA:18691","RHEA:18692")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPC", replacement = df_reaction[["sn2LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## pc_to_pa
    if (reaction %in% c("RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PA", replacement = df_reaction[["PA"]])
    }
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PS", replacement = df_reaction[["PS"]])
    }
    
    ## pco_to_lpco
    if (reaction %in% c("RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC-O", replacement = df_reaction[["PCO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPC-O", replacement = df_reaction[["LPCO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    } 
    
    if (reaction == "pe_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## pe_to_sn1lpe
    if (reaction %in% c("RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPE", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## pe_to_sn2lpe
    if (reaction %in% c("RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPE", replacement = df_reaction[["sn2LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "pe_to_nape_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPE", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPC", replacement = df_reaction[["sn2LPC"]])
    }
    
    if (reaction == "pe_to_nape_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPE", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])   
    }
    
    if (reaction == "pe_to_pa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PA", replacement = df_reaction[["PA"]])
    }

    ## pe_to_ps
    if (reaction %in% c("RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PS", replacement = df_reaction[["PS"]])
    }
    
    if (reaction == "peo_to_lpeo") {
        # adjust variable substrates and products
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPE-O", replacement = df_reaction[["LPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "peo_to_napeo_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPEO", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPC", replacement = df_reaction[["sn2LPC"]])
    }
    
    if (reaction == "peo_to_napeo_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPEO", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])
    }
    
    ## peo_to_pep
    if (reaction %in% c("RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-O", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-P", replacement = df_reaction[["PEP"]])
    }
    
    ## pep_to_lpep
    if (reaction %in% c("RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-P", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_LPE-P", replacement = df_reaction[["LPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "pep_to_napep_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-P", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPEP", replacement = df_reaction[["NAPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2-LPC", replacement = df_reaction[["sn2LPC"]])
    }
    
    if (reaction == "pep_to_napep_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PC", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE-P", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_NAPEP", replacement = df_reaction[["NAPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1-LPC", replacement = df_reaction[["sn1LPC"]])
    }
    
    ## pg_to_cl
    if (reaction %in% c("RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CDP-DG", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PG", replacement = df_reaction[["PG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CL", replacement = df_reaction[["CL"]])
    }
    
    ## pgp_to_pg
    if (reaction %in% c("RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pgp", replacement = df_reaction[["PGP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PG", replacement = df_reaction[["PG"]])
    }
    
    ## pi_to_dg
    if (reaction %in% c("RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PI", replacement = df_reaction[["PI"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1,2-DG", replacement = df_reaction[["DG"]])
    }
    
    ## pi_to_sn1lpi
    if (reaction %in% c("RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_PI", replacement = df_reaction[["PI"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1-LPI", replacement = df_reaction[["sn1LPI"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_FA", replacement = df_reaction[["FA"]])
    }
    
    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PS", replacement = df_reaction[["PS"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_PE", replacement = df_reaction[["PE"]])
    }
    
    if (reaction == "sm_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_CER", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_SM", replacement = df_reaction[["SM"]])
    }
    
    if (reaction == "sphinga_to_dhcer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_AcylCoA", replacement = df_reaction[["AcylCoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_DHCER", replacement = df_reaction[["DHCER"]])
    }
    
    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274",
            "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867")) {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_TG",
            replacement = c(df_reaction[["TG"]], df_reaction[["TG"]]))
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1,2-DG", replacement = c(df_reaction[["sn1Loss_DG"]], 
                df_reaction[["sn3Loss_DG"]]))
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_FA", replacement = c(df_reaction[["sn1Loss_FA"]], 
                df_reaction[["sn3Loss_FA"]]))
    }
    
    ## adjust fixed substrates and products
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CoA ", replacement = "CoA ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CoA$", replacement = "CoA")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CoA_c ", replacement = "CoA_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CoA_c$", replacement = "CoA_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CMP ", replacement = "CMP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CMP$", replacement = "CMP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CTP ", replacement = "CTP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CTP$", replacement = "CTP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_AMP ", replacement = "AMP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_AMP$", replacement = "AMP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ADP ", replacement = "ADP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ADP$", replacement = "ADP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ATP ", replacement = "ATP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ATP$", replacement = "ATP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ADP ", replacement = "ADP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ADP$", replacement = "ADP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_UDP ", replacement = "UDP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_UDP$", replacement = "UDP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_UDP-Glucose ", replacement = "UDP-Glucose ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_UDP-Glucose$", replacement = "UDP-Glucose")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_NADH ", replacement = "NADH ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_NADH$", replacement = "NADH")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_NAD ", replacement = "NAD+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_NAD$", replacement = "NAD+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_NADPH ", replacement = "NADPH ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_NADPH$", replacement = "NADPH")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_NADP ", replacement = "NADP+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_NADP$", replacement = "NADP+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_O2 ", replacement = "O2 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_O2$", replacement = "O2")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H2O ", replacement = "H2O ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H2O_c$", replacement = "H2O_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H2O_c ", replacement = "H2O_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H2O$", replacement = "H2O")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CO2 ", replacement = "CO2 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CO2$", replacement = "CO2")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H+ ", replacement = "H+ ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H+_c ", replacement = "H+_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H+$", replacement = "H+")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H+_c$", replacement = "H+_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H+ ", replacement =  "H+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H+$", replacement =  "H+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Pi ", replacement = "Pi ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Pi$", replacement = "Pi")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_PPi ", replacement = "PPi ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_PPi$", replacement = "PPi")
    .formula <- stringi::stri_replace_all_fixed(str = .formula,
        pattern = "M_Fe2+-cytochrome_b5 ", replacement = "Fe2+-cytochrome b5 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula,
        pattern = "M_Fe2+-cytochrome_b5$", replacement = "Fe2+-cytochrome b5$")
    .formula <- stringi::stri_replace_all_fixed(str = .formula,
        pattern = "M_Fe3+-cytochrome_b5 ", replacement = "Fe3+-cytochrome b5 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula,
        pattern = "M_Fe3+-cytochrome_b5$", replacement = "Fe3+-cytochrome b5$")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Glycerol ", replacement = "Glycerol ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Glycerol$", replacement = "Glycerol")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Glycerol-3-P ", replacement = "Glycerol-3-P ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Glycerol-3-P$", replacement = "Glycerol-3-P")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_myo-Inositol ", replacement = "myo-Inositol ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_myo-Inositol$", replacement = "myo-Inositol")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_myo-Inositol-1-P ", replacement = "myo-Inositol-1-P ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_myo-Inositol-1-P$", replacement = "myo-Inositol-1-P")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Dihydroxyacetone-P ", replacement = "Dihydroxyacetone-P ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Dihydroxyacetone-P$", replacement = "Dihydroxyacetone-P")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CDP-Choline ", replacement = "CDP-Choline ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CDP-Choline$", replacement = "CDP-Choline")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_L-Serine_c ", replacement = "L-Serine_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_L-Serine_c$", replacement = "L-Serine_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_L-Serine ", replacement = "L-Serine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_L-Serine$", replacement = "L-Serine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Ethanolamine ", replacement = "Ethanolamine ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Ethanolamine_c ", replacement = "Ethanolamine_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Ethanolamine$", replacement = "Ethanolamine")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Ethanolamine_c$", replacement = "Ethanolamine_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_P-Ethanolamine ", replacement = "P-Ethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_P-Ethanolamine$", replacement = "P-Ethanolamine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_CDP-Ethanolamine ", replacement = "CDP-Ethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_CDP-Ethanolamine$", replacement = "CDP-Ethanolamine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_c17isosphgn ", replacement = "SPH(d16:0(1OH,3OH)(15Me)) ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_c17isosphgn$", replacement = "SPH(d16:0(1OH,3OH)(15Me))")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Choline_c ", replacement = "Choline_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Choline_c$", replacement = "Choline_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula,
        pattern = "M_Choline ", replacement = "Choline ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Choline$", replacement = "Choline")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_P-Choline ", replacement = "P-Choline ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_P-Choline$", replacement = "P-Choline")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Glycerophosphoethanolamine ", replacement = "Glycerophosphoethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Glycerophosphoethanolamine$", replacement = "Glycerophosphoethanolamine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_Glycerophosphocholine ", replacement = "Glycerophosphocholine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_Glycerophosphocholineg3pc$", replacement = "Glycerophosphocholine")
    
    ## replace existing lipid backbones (DG, PC)
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_1,2-DG ", replacement = "DG ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_1,2-DG$", replacement = "DG")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_PC ", replacement = "PC ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_PC$", replacement = "PC")

    ## write back to entry reaction_formula
    df$reaction_formula <- .formula
    
    ## split the df$reaction_formula into substrates and products 
    .formula_tmp <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = ">|<", replacement = "")
    .formula_tmp <- stringi::stri_split(str = .formula_tmp, fixed = " = ")
    .formula_substrate <- unlist(lapply(.formula_tmp, "[", 1))
    .formula_product <- unlist(lapply(.formula_tmp, "[", 2))

    ## write the substrates and products to the df
    df$reaction_substrate <- .formula_substrate
    df$reaction_product <- .formula_product
    
    ## return the data.frame
    df
}