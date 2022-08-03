
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
    
    if (reaction %in% "tg_to_dg") .char <- c(.char, .char)
    
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
    if (reaction == "RHEA:36171") { ## acdhap_to_alkyldhap
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_adhap", replacement = df_reaction[["ACDHAP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkylR1oh", replacement = df_reaction[["FATOH"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akdhap", replacement <- df_reaction[["ALKYLDHAP"]])
    }
    
    if (reaction == "RHEA:36175") { ## alkyldhap_to_lpao
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akdhap", replacement = df_reaction[["ALKYLDHAP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alpa_pl", replacement = df_reaction[["LPAO"]])
    }
    
    if (reaction == "c1p_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrmp", replacement = df_reaction[["C1P"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
    }
    
    if (reaction == "RHEA:12593") { ## cdpdg_to_pgp
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_cdpdag", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pgp", replacement = df_reaction[["PGP"]])
    }
    
    if (reaction == "RHEA:11580") { ## cdpdg_to_pi
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_cdpdag", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pail", replacement = df_reaction[["PI"]])
    }
    
    if (reaction == "RHEA:17929") { ## cer_to_c1p
        .formula <- stringi::stri_replace_first_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrmp", replacement = df_reaction[["C1P"]])
    }
    
    if (reaction == "RHEA:12088") { ## cer_to_glccer
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isogluside", replacement = df_reaction[["GLCCER"]])
    }
    
    if (reaction == "cer_to_sm") { ## 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isosphmyln", replacement = df_reaction[["SM"]])
    }
    
    if (reaction == "coa_to_acdhap") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_adhap", replacement = df_reaction[["DHAP"]])
    }
    
    if (reaction == "RHEA:52716") { ## coa_to_fatoh
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkylR1oh",  replacement = df_reaction[["FATOH"]])
    }
    
    if (reaction == "RHEA:15325") { ## coa_to_lpa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alpa_pl", replacement = df_reaction[["LPA"]])
    }
    
    if (reaction == "dg_to_sn1mg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
                pattern = "M_12dag", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_1magol", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "dg_to_sn2mg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_mag", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:10272") { ## dg_to_pa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "RHEA:32939") { ## dg_to_pc
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
    }
    
    if (reaction == "RHEA:32943") { ## dg_to_pe
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
    }
    
    if (reaction == "RHEA:10868") { ## dg_to_tg
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_tag", replacement = df_reaction[["TG"]])
    }
    
    if (reaction == "RHEA:36179") { ## dgo_to_pco
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_cdpchol", replacement = "CDP-Chol")
        .formula <- stringi::stri_replace_first_fixed(str = .formula, 
            pattern = "M_akac2g", replacement = df_reaction[["DGO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gchol", replacement = df_reaction[["PCO"]])
    }

    if (reaction == "RHEA:36187") { ## dgo_to_peo
        .formula <- stringi::stri_replace_first_fixed(str = .formula,
            pattern = "M_akac2g", replacement = df_reaction[["DGO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
    }
    
    if (reaction == "dhcer_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhcrm", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
    }
    
    if (reaction == "dhcer_to_dhsm") { ## 
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhcrm", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhsphmyln", replacement = df_reaction[["DHSM"]])
    }
    
    if (reaction == "RHEA:19253") { ## dhsm_to_dhcer
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhcrm", replacement = df_reaction[["DHCER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhsphmyln", replacement = df_reaction[["DHSM"]])
    }
    
    if (reaction == "RHEA:15421") { ## fa_to_coa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
    }
    
    if (reaction == "RHEA:45420") { ## lnape_to_gpnae
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_lnape", replacement = df_reaction[["LNAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_gpnae", replacement = df_reaction[["GPNAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:19709") { ## lpa_to_pa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alpa_pl", replacement = df_reaction[["LPA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
    }
    
    ## if (reaction == "RHEA:19709") { ## lpao_to_pao
    ##     .formula <- stringi::stri_replace_all_fixed(str = .formula,
    ##         pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
    ##     .formula <- stringi::stri_replace_all_fixed(str = .formula, 
    ##         pattern = "M_alpa_pl", replacement = df_reaction[["LPAO"]])
    ##     .formula <- stringi::stri_replace_all_fixed(str = .formula, 
    ##         pattern = "M_pa_pl", replacement = df_reaction[["PAO"]])
    ## }
    
    if (reaction == "RHEA:15177") { ## sn1lpc_to_fa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pc", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:44696") { ## sn2lpc_to_fa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpc", replacement = df_reaction[["sn2LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:12937") { ## sn1lpc_to_pc
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pc", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
    }
    
    if (reaction == "RHEA:32967") { ## sn1lpe_to_fa
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_ag3pe", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "sn2lpe_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpe", replacement = df_reaction[["sn2LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:32995") { ## sn1lpe_to_pe
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_acg3pe", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
    }
    
    if (reaction == "lpeo_to_peo") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ak2lgpe", replacement = df_reaction[["LPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
    }
    
    if (reaction == "lpep_to_pep") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alken2gpe", replacement = df_reaction[["LPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2gpe", replacement = df_reaction[["PEP"]])
    }
    
    if (reaction == "sn1mg_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1magol", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "sn2mg_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_mag", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "sn1mg_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1magol", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])   
    }
    
    if (reaction == "sn2mg_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_mag", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "sn1mg_to_lpa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1magol", replacement = df_reaction[["sn1MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alpa_pl", replacement = df_reaction[["LPA"]])
    }
    
    if (reaction == "sn2mg_to_sn1mg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_mag", replacement = df_reaction[["sn2MG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_1magol", replacement = df_reaction[["sn1MG"]])
    }
    
    if (reaction == "nae_to_fa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_nae", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "nape_to_lnape") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_nape", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_lnape", replacement = df_reaction[["LNAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "nape_to_nae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_nape", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_nae", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
    }
    
    if (reaction == "nape_to_pnae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_nape", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_pnae", replacement = df_reaction[["PNAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "napeo_to_nae") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_akac2nape", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_nae", replacement = df_reaction[["NAE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_akac2gp", replacement = df_reaction[["PAO"]])
    }
    
    if (reaction == "pa_to_cdpdg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_cdpdag", replacement = df_reaction[["CDPDG"]])
    }
    
    if (reaction == "pa_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "pao_to_dgo") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_akac2gp", replacement = df_reaction[["PAO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2g", replacement = df_reaction[["DGO"]])
    }
    
    if (reaction == "RHEA:10604") { ## pc_to_dg
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "RHEA:15801") { ## pc_to_sn1lpc
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula,
            pattern = "M_ag3pc", replacement = df_reaction[["sn1LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:18689") { ## pc_to_sn2lpc
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpc", replacement = df_reaction[["sn2LPC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:14445") { ## pc_to_pa
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
    }
    
    if (reaction == "RHEA:45088") { ## pc_to_ps
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pchol", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ps", replacement = df_reaction[["PS"]])
    }
    
    if (reaction == "pco_to_lpco"){
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gchol", replacement = df_reaction[["PCO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ak2lgchol", replacement = df_reaction[["LPCO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    } 
    
    if (reaction == "pe_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = df_reaction[["DG"]])
    }
    
    if (reaction == "pe_to_sn1lpe") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pe", replacement = df_reaction[["sn1LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "RHEA:44408") { ## pe_to_sn2lpe
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpe", replacement = df_reaction[["sn2LPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "pe_to_nape_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_nape", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpc", replacement = df_reaction[["LPC"]])
    }
    
    if (reaction == "pe_to_nape_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_nape", replacement = df_reaction[["NAPE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pc", replacement = df_reaction[["LPC"]])   
    }
    
    if (reaction == "pe_to_pa") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pa_pl", replacement = df_reaction[["PA"]])
    }
    
    if (reaction == "RHEA:27606") { ## pe_to_ps
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ps", replacement = df_reaction[["PS"]])
    }
    
    if (reaction == "peo_to_lpeo") {
        # adjust variable substrates and products
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ak2lgpe", replacement = df_reaction[["LPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "peo_to_napeo_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2nape", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpc", replacement = df_reaction[["LPC"]])
    }
    
    if (reaction == "peo_to_napeo_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2nape", replacement = df_reaction[["NAPEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pc", replacement = df_reaction[["LPC"]])
    }
    
    if (reaction == "peo_to_pep") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_akac2gpe", replacement = df_reaction[["PEO"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2gpe", replacement = df_reaction[["PEP"]])
    }
    
    if (reaction == "pep_to_lpep") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2gpe", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alken2gpe", replacement = df_reaction[["LPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = df_reaction[["FA"]])
    }
    
    if (reaction == "pep_to_napep_sn1") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2gpe", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2nape", replacement = df_reaction[["NAPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_2agpc", replacement = df_reaction[["LPC"]])
    }
    
    if (reaction == "pep_to_napep_sn2") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pc", replacement = df_reaction[["PC"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2gpe", replacement = df_reaction[["PEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_alkenac2nape", replacement = df_reaction[["NAPEP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ag3pc", replacement = df_reaction[["LPC"]])
    }
    
    if (reaction == "pg_to_cl") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_cdpdag", replacement = df_reaction[["CDPDG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pg", replacement = df_reaction[["PG"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_clpn", replacement = df_reaction[["CL"]])
    }
    
    if (reaction == "pgp_to_pg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pgp", replacement = df_reaction[["PGP"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pg", replacement = df_reaction[["PG"]])
    }
    
    if (reaction == "ps_to_pe") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_ps", replacement = df_reaction[["PS"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_pe", replacement = df_reaction[["PE"]])
    }
    
    if (reaction == "sm_to_cer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isocrm", replacement = df_reaction[["CER"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isosphmyln", replacement = df_reaction[["SM"]])
    }
    
    if (reaction == "sphinga_to_dhcer") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fataccoa", replacement = df_reaction[["CoA"]])
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_c17isodhcrm", replacement = df_reaction[["DHCER"]])
    }
    
    if (reaction == "tg_to_dg") {
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_tag",
            replacement = c(df_reaction[["TG"]], df_reaction[["TG"]]))
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_12dag", replacement = c(df_reaction[["sn1Loss_dg"]], 
                df_reaction[["sn3Loss_dg"]]))
        .formula <- stringi::stri_replace_all_fixed(str = .formula, 
            pattern = "M_fatacid", replacement = c(df_reaction[["sn1Loss_fa"]], 
                df_reaction[["sn3Loss_fa"]]))
    }
    
    ## adjust fixed substrates and products
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_coa ", replacement = "CoA ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_coa$", replacement = "CoA")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_cmp ", replacement = "CMP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_cmp$", replacement = "CMP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ctp ", replacement = "CTP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ctp$", replacement = "CTP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_amp ", replacement = "AMP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_amp$", replacement = "AMP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_adp ", replacement = "ADP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_adp$", replacement = "ADP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_atp ", replacement = "ATP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_atp$", replacement = "ATP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_adp ", replacement = "ADP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_adp$", replacement = "ADP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_udp ", replacement = "UDP ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_udp$", replacement = "UDP")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_udpg ", replacement = "UDP-Glucose ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_udpg$", replacement = "UDP-Glucose")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_nadh ", replacement = "NADH ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_nadh$", replacement = "NADH")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_nad ", replacement = "NAD+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_nad$", replacement = "NAD+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_nadph ", replacement = "NADPH ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_nadph$", replacement = "NADPH")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_nadp ", replacement = "NADP+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_nadp$", replacement = "NADP+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_o2 ", replacement = "O2 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_o2$", replacement = "O2")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_h2o ", replacement = "H2O ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_h2o_c$", replacement = "H2O_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_h2o_c ", replacement = "H2O_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_h2o$", replacement = "H2O")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_co2 ", replacement = "CO2 ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_co2$", replacement = "CO2")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_h ", replacement = "H+ ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_h_c ", replacement = "H+_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_h$", replacement = "H+")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_h_c$", replacement = "H+_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_H ", replacement =  "H+ ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_H$", replacement =  "H+")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_pi ", replacement = "Pi ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_pi$", replacement = "Pi")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ppi ", replacement = "PPi ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ppi$", replacement = "PPi")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_glyc ", replacement = "Glycerol ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_glyc$", replacement = "Glycerol")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_glyc3p ", replacement = "Glycerol-3-P ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_glyc3p$", replacement = "Glycerol-3-P")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_inost ", replacement = "Inositol ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_inost$", replacement = "Inositol")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_pchol ", replacement = "PC ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_pchol$", replacement = "PC")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_12dag ", replacement = "DG ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_12dag$", replacement = "DG")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_dhap ", replacement = "Dihydroxyacetone-P ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_dhap$", replacement = "Dihydroxyacetone-P")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_cdpea ", replacement = "CDP-Ethn ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_cdpea$", replacement = "CDP-Ethn")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_cdpchol ", replacement = "CDP-Chol ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_cdpchol$", replacement = "CDP-Chol")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_12dag ", replacement = "DG ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_12dag$", replacement = "DG")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ser_L_c ", replacement = "L-Serine_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ser_L_c$", replacement = "L-Serine_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ser_L ", replacement = "L-Serine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ser_L$", replacement = "L-Serine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_etha ", replacement = "Ethanolamine ")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_etha_c ", replacement = "Ethanolamine_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_etha$", replacement = "Ethanolamine")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_etha_c$", replacement = "Ethanolamine_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_ethamp ", replacement = "P-Ethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_ethamp$", replacement = "P-Ethanolamine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_c17isosphgn ", replacement = "SPH(d16:0(1OH,3OH)(15Me)) ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_c17isosphgn$", replacement = "SPH(d16:0(1OH,3OH)(15Me))")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_chol_c ", replacement = "Choline_c ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_chol_c$", replacement = "Choline_c")
    .formula <- stringi::stri_replace_all_fixed(str = .formula,
        pattern = "M_chol ", replacement = "Choline ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_chol$", replacement = "Choline")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_cholp ", replacement = "P-Choline ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_cholp$", replacement = "P-Choline")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_g3pe ", replacement = "Glycerophosphoethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_g3pe$", replacement = "Glycerophosphoethanolamine")
    .formula <- stringi::stri_replace_all_fixed(str = .formula, 
        pattern = "M_g3pc ", replacement = "Glycerophosphoethanolamine ")
    .formula <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = "M_g3pc$", replacement = "Glycerophosphoethanolamine")

    ## write back to entry reaction_formula
    df$reaction_formula <- .formula
    
    ## split the df$reaction_formula into substrates and products 
    .formula_tmp <- stringi::stri_replace_all_regex(str = .formula, 
        pattern = ">|<", replacement = "")
    .formula_tmp <- stringi::stri_split(str = .formula_tmp, fixed = "=")
    .formula_substrate <- unlist(lapply(.formula_tmp, "[", 1))
    .formula_product <- unlist(lapply(.formula_tmp, "[", 2))
    
    ## remove the trailing spaces for substrates and products
    #.formula_substrate <- stringi::stri_replace_all_regex(str = .formula_substrate, 
    #    pattern = "[+]", replacement = ",")
    .formula_substrate <- stringi::stri_replace_all_fixed(str = .formula_substrate, 
        pattern = " ", replacement = "")
    #.formula_product <- stringi::stri_replace_all_regex(str = .formula_product, 
    #    pattern = "[+]", replacement = ",")
    .formula_product <- stringi::stri_replace_all_fixed(str = .formula_product, 
        pattern = " ", replacement = "")
    
    ## write the substrates and products to the df
    df$reaction_substrate <- .formula_substrate
    df$reaction_product <- .formula_product
    
    ## return the data.frame
    df
}