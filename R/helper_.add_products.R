#' @name .add_products
#' 
#' @title Produce products by replacing strings of substrates 
#' 
#' @description 
#' Helper function for \code{create_reaction}.
#' 
#' The function \code{.add_products} returns a data.frame of substrates and
#' products. The products are taken and modified from the substrates.
#' 
#' @details
#' The string replacement depend on the \code{reaction} argument.
#' 
#' Depending on the \code{argument}, different columns (i.e. substrates)
#' are required and different columns (i.e. products) are returned.
#' 
#' The function calls internally the function
#' \code{stringi::stri_replace_all_fixed}.
#' 
#' @param substrate data.frame
#' @param reaction \code{character(1)}
#' 
#' @return data.frame
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'    and Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' 
#' ## create template
#' template <- wormLipidPredictR:::.create_template(template = NA, 
#'     reaction = "fa_to_coa")
#' 
#' ## create data.frame of substrates
#' df_substrates <- wormLipidPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#'     
#' ## add products to data.frame
#' wormLipidPredictR:::.add_products(substrates = df_substrates, 
#'     reaction = "fa_to_coa")
#'  
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
.add_products <- function(substrates, reaction = "fa_to_coa") {

    ## make the substrates argument a little bit smaller
    .s <- substrates
    ## start of if statements 
    ## depending on the reaction type, replace the expressions of the substrates
    ## and products in the reaction
    if (reaction == "acdhap_to_alkyldhap") {
        .s$FA <- stringi::stri_replace_all_regex(str = .s$ACDHAP, 
            pattern = "DHAP\\(", replacement = "FA(")
        .s$ALKYLDHAP <- stringi::stri_replace_all_regex(str = .s$FATOH, 
            pattern = "FATOH\\(", replacement = "DHAP(O-")
    }

    if (reaction == "alkyldhap_to_lpao") {
        .s$LPAO <- stringi::stri_replace_all_fixed(str = .s$ALKYLDHAP, 
            pattern = "DHAP", replacement = "PA")
        .s$LPAO <- stringi::stri_replace_all_regex(str = .s$LPAO, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }

    if (reaction == "c1p_to_cer") {
        .s$CER <- stringi::stri_replace_all_fixed(str = .s$C1P, 
            pattern = "C1P", replacement = "Cer")
    }
    
    if (reaction == "cdpdg_to_pgp") {
        .s$PGP <- stringi::stri_replace_all_fixed(str = .s$CDPDG, 
            pattern = "CDP-DG", replacement = "PGP")
    }
            
    if (reaction == "cdpdg_to_pi") {
        .s$PI <- stringi::stri_replace_all_fixed(str = .s$CDPDG, 
            pattern = "CDP-DG", replacement = "PI")       
    }
            
    if (reaction == "cer_to_c1p") {
        .s$C1P <- stringi::stri_replace_all_fixed(str = .s$CER, 
            pattern = "Cer", replacement = "C1P")
    }
            
    if (reaction == "cer_to_glccer") {
        .s$GLCCER <- stringi::stri_replace_all_fixed(str = .s$CER, 
            pattern = "Cer", replacement = "GlcCer")
    }
            
    if (reaction == "cer_to_sm") {
        .s$SM <- stringi::stri_replace_all_fixed(str = .s$CER, 
            pattern = "Cer", replacement = "SM")
    }
    
    if (reaction == "coa_to_acdhap") {
        .s$DHAP <- stringi::stri_replace_all_fixed(str = .s$CoA, 
            pattern = "CoA", replacement = "DHAP")
    }
    
    if (reaction == "coa_to_fatoh") {
        .s$FATOH <- stringi::stri_replace_all_fixed(str = .s$CoA, 
            pattern = "CoA", replacement = "FATOH")
    }
    
    if (reaction == "coa_to_lpa") {
        .s$LPA <- stringi::stri_replace_all_fixed(str = .s$CoA, pattern = "CoA", 
            replacement = "PA")
        .s$LPA <- stringi::stri_replace_all_regex(str = .s$LPA, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }

    if (reaction == "dg_to_sn1mg") {
        .s$sn1MG <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$DG), 
            function(f) {paste0("MG(", f[1], "/0:0/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$DG), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    if (reaction == "dg_to_sn2mg") {
        .s$sn2MG <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$DG), 
            function(f) {paste0("MG(0:0/", f[2], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$DG), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
            
    if (reaction == "dg_to_pa") {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PA")
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$PA, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    if (reaction == "dg_to_pc") {
        .s$PC <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PC")
        .s$PC <- stringi::stri_replace_all_fixed(str = .s$PC, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    if (reaction == "dg_to_pe") {
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PE")
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$PE, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }
            
    if (reaction == "dg_to_tg") {
        .s$TG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "/0:0", 
            replacement = paste0("/", 
                unlist(lipidomicsUtils::isolate_radyls(.s$CoA))))
         .s$TG <- stringi::stri_replace_all_fixed(str = .s$TG, pattern = "^DG", 
            replacement = "TG")
    }
    
    if (reaction == "dgo_to_pco") {
        .s$PCO <- stringi::stri_replace_all_fixed(str = .s$DGO, pattern = "DG", 
            replacement = "PC")
        .s$PCO <- stringi::stri_replace_all_regex(str = .s$PCO, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    if (reaction == "dgo_to_peo") {
        .s$PEO <- stringi::stri_replace_all_fixed(str = .s$DGO, pattern = "DG", 
            replacement = "PE")
        .s$PEO <- stringi::stri_replace_all_fixed(str = .s$PEO, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }
            
    if (reaction == "dhcer_to_cer") {
        .s$CER <- stringi::stri_replace_all_fixed(str = .s$DHCER, 
            pattern = "Cer\\(d16:0\\(3OH,4OH\\)\\(15Me\\)\\/",
            replacement = "Cer(d16:1(4E)(3OH,4OH)(15Me)/")
    }
            
    if (reaction == "dhcer_to_dhsm") {
        .s$DHSM <- stringi::stri_replace_all_fixed(str = .s$DHCER, 
            pattern = "Cer", replacement = "SM")
    }
            
        
    if (reaction == "dhsm_to_dhcer") {
        .s$DHCER <- stringi::stri_replace_all_fixed(str = .s$DHSM, 
            pattern = "SM", replacement = "Cer")
    }
        
    if (reaction == "fa_to_coa") {
        .s$CoA <- stringi::stri_replace_all_fixed(str = .s$FA, "FA", 
            replacement = "CoA")
    }
    
    if (reaction == "lnape_to_gpnae") {
        .s$GPNAE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$LNAPE), 
            function(f) {paste0("GPNAE(", f[3], ")")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$LNAPE), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
        
    if (reaction == "lpa_to_pa") {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$LPA, pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f)})))
    }
            
    if (reaction == "lpao_to_pao") {
        .s$PAO <- stringi::stri_replace_all_fixed(str = .s$LPAO, 
            pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f)})))
    }

    if (reaction == "sn1lpc_to_fa") {
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn1LPC), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    if (reaction == "sn2lpc_to_fa") {
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn2LPC), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "sn1lpc_to_pc") {
        .s$PC <- stringi::stri_replace_all_fixed(str = .s$sn1LPC, 
            pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f)})))
    }

    if (reaction == "sn1lpe_to_fa") {
        ## sn2 loss
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn1LPE), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    if (reaction == "sn2lpe_to_fa") {
        ## sn2 loss 
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn2LPE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "sn1lpe_to_pe") {
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$sn1LPE, 
            pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA),
                function(f) {paste0("/", f)})))
    }

    if (reaction == "lpeo_to_peo") {
        .s$PEO <- stringi::stri_replace_all_fixed(str = .s$LPEO, 
            pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f)})))
    }

    if (reaction == "lpep_to_pep") {
        .s$PEP <- stringi::stri_replace_all_fixed(str = .s$LPEP, 
            pattern = "/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f)})))
    }

    if (reaction == "sn1mg_to_dg") {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$sn1MG, 
            pattern = "/0:0/0:0", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("/", f, "/0:0")})))
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "MG",
            replacement = "DG")
    }
    
    if (reaction == "sn2mg_to_dg") {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$sn2MG, 
            pattern = "\\(0:0/", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$CoA), 
                function(f) {paste0("(", f, "/")})))
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "MG", 
            replacement = "DG")
    }
    
    if (reaction == "sn1mg_to_fa") {
        ## "sn2 loss"
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn1MG), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    if (reaction == "sn2mg_to_fa") {
        ## "sn2 loss"
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn2MG), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "sn1mg_to_lpa") {
        .s$LPA <- stringi::stri_replace_all_fixed(str = .s$sn1MG, 
            pattern = "MG", replacement = "PA")
        .s$LPA <- stringi::stri_replace_all_fixed(str = .s$LPA, "/0:0\\)$", 
            replacement = "\\)")
    }

    if (reaction == "sn2mg_to_sn1mg") {
        ## "sn2 loss"
        .s$sn1MG <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$sn2MG), 
            function(f) {paste0("MG(", f[2], "/0:0/0:0)")}))
    }

    if (reaction == "nae_to_fa") {
        .s$FA <- stringi::stri_replace_all_fixed(str = .s$NAE, pattern = "NAE", 
            replacement = "FA")
    }

    if (reaction == "nape_to_lnape") {
        .s$LNAPE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("NAPE(", f[1], "/0:0/", f[3], ")")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "nape_to_nae") {
        .s$NAE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("NAE(", f[3], ")")}))
        .s$PA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("PA(", f[1], "/", f[2], ")")}))
    }

    if (reaction == "nape_to_pnae") {
        .s$PNAE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("PNAE(", f[3], ")")}))
        .s$DG <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPE), 
            function(f) {paste0("DG(", f[1], "/", f[2], "/0:0)")}))
    }
    
    if (reaction == "napeo_to_nae") {
        .s$NAE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPEO), 
            function(f) {paste0("NAE(", f[3], ")")}))
        .s$PAO <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$NAPEO), 
            function(f) {paste0("PA(", f[1], "/", f[2], ")")}))
    }
    
    if (reaction == "pa_to_cdpdg") {
        .s$CDPDG <- stringi::stri_replace_all_fixed(str = .s$PA, pattern = "PA",
            replacement = "CDP-DG")
    }

    if (reaction == "pa_to_dg") {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PA, pattern = "PA", 
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "pao_to_dgo") {
        .s$DGO <- stringi::stri_replace_all_fixed(str = .s$PAO, pattern = "PA", 
            replacement = "DG")
        .s$DGO <- stringi::stri_replace_all_regex(str = .s$DGO, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }

    if (reaction == "pc_to_dg") {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC", 
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "pc_to_sn1lpc") {
        ## "sn2 loss"
        .s$sn1LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "pc_to_sn2lpc") {
        ## "sn2 loss"
        .s$sn2LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
    
    if (reaction == "pc_to_pa") {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC", 
            replacement = "PA")
    }
    
    if (reaction == "pc_to_ps") {
        .s$PS <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC",
            replacement = "PS")
    }
    
    if (reaction == "pco_to_lpco") {
        ## "sn2 loss"
        .s$LPCO <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PCO), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PCO), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    if (reaction == "pe_to_dg") {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE",
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$",
            replacement = "/0:0\\)")
    }
    
    if (reaction == "pe_to_sn1lpe") {
        ## "sn2" loss
        .s$sn1LPE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PE), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "pe_to_sn2lpe") {
        .s$sn2LPE <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PE), 
            function(f) {paste0("PE(0:0/", f[1], ")")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    if (reaction == "pe_to_nape_sn1") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPE <- stringi::stri_replace_all_regex(str = .s$PE, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPE <- stringi::stri_replace_all_fixed(str = .s$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }

    if (reaction == "pe_to_nape_sn2") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPE <- stringi::stri_replace_all_regex(str = .s$PE, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPE <- stringi::stri_replace_all_fixed(str = .s$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }
    
    if (reaction == "pe_to_pa") {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE", 
            replacement = "PA")
    }
    
    if (reaction == "pe_to_ps") {
        .s$PS <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE", 
            replacement = "PS")
    }
    
    if (reaction == "peo_to_lpeo") {
        ## "sn2 loss"
        .s$LPEO <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PEO), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PEO), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "peo_to_napeo_sn1") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPEO <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPEO <- stringi::stri_replace_all_fixed(str = .s$NAPEO, 
            pattern = "PE", replacement = "NAPE")
    }
    
    if (reaction == "peo_to_napeo_sn2") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPEO <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPEO <- stringi::stri_replace_all_fixed(str = .s$NAPEO, 
            pattern = "PE", replacement = "NAPE")
    }
    
    if (reaction == "peo_to_pep") {
        .s$PEP <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "PE\\(O-", replacement = "PE(P-")
    }
    
    if (reaction == "pep_to_lpep") {
        ## "sn2 loss"
        .s$LPEP <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PEP), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PEP), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    if (reaction == "pep_to_napep_sn1") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPEP <- stringi::stri_replace_all_regex(str = .s$PEP, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPEP <- stringi::stri_replace_all_fixed(str = .s$NAPEP,
            pattern = "PE", replacement = "NAPE")
    }

    if (reaction == "pep_to_napep_sn2") {
        .s$LPC <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPEP <- stringi::stri_replace_all_regex(str = .s$PEP, 
            pattern = "\\)$", 
            replacement = unlist(lapply(lipidomicsUtils::isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPEP <- stringi::stri_replace_all_fixed(str = .s$NAPEP, 
            pattern = "PE", replacement = "NAPE")
    }
    
    if (reaction == "pg_to_cl") {
        ## isolate core from PGs
        .s$PGs2 <- stringi::stri_replace_all_fixed(str = .s$PG, 
            pattern = "^PG\\(", replacement = "")
        .s$PGs2 <- stringi::stri_replace_all_regex(str = .s$PGs2, 
            pattern = "\\)$", replacement = "")
         
        ## isolate core from CDPDGs
        .s$CDPDGs2 <- stringi::stri_replace_all_fixed(str = .s$CDPDG, 
            pattern = "^CDP-DG\\(", replacement = "")
        .s$CDPDGs2 <- stringi::stri_replace_all_regex(str = .s$CDPDGs2, 
            pattern = "\\)$", replacement = "")
        
        ## create CL
        .s$CL <- paste0("CL(1'-[", .s$PGs2, "],3'-[", .s$CDPDGs2, "])")
    }
    
    if (reaction == "pgp_to_pg") {
        .s$PG <- stringi::stri_replace_all_fixed(str = .s$PGP, pattern = "PGP", 
            replacement = "PG")
    }

    if (reaction == "ps_to_pe") {
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$PS, pattern = "PS",
            replacement = "PE")
    }
    
    if (reaction == "sm_to_cer") {
        .s$CER <- stringi::stri_replace_all_fixed(str = .s$SM, pattern = "SM",
            replacement = "Cer")
    }

    if (reaction == "sphinga_to_dhcer") {
        .s$DHCER <- stringi::stri_replace_all_regex(str = .s$CoA, 
            pattern = "CoA\\(", replacement = "Cer(d16:0(3OH,4OH)(15Me)/")
    }

    if (reaction == "tg_to_dg") {
        ## "sn1 loss"
        .s$sn1Loss_dg <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$TG), 
            function(f) {paste0("DG(", f[3], "/", f[2], "/0:0)")}))
        .s$sn1Loss_fa <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$TG), 
            function(f) {paste0("FA(", f[1], ")")}))

        ## "sn3 loss"
        .s$sn3Loss_dg <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$TG), 
            function(f) {paste0("DG(", f[1], "/", f[2], "/0:0)")}))
        .s$sn3Loss_fa <- unlist(lapply(lipidomicsUtils::isolate_radyls(.s$TG), 
            function(f) {paste0("FA(", f[3], ")")}))
    }
    
    ## end of if statements
    ## return the data frame with added product(s)
    .s
}


