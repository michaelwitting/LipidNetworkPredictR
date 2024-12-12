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
#' @param substrates \code{data.frame}
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
#' template <- LipidNetworkPredictR:::.create_template(template = list(), 
#'     reaction = "RHEA:15421")
#' 
#' ## create data.frame of substrates
#' df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     template = template)
#'     
#' ## add products to data.frame
#' LipidNetworkPredictR:::.add_products(substrates = df_substrates, 
#'     reaction = "RHEA:15421")
#'  
#' @importFrom stringi stri_replace_all_fixed stri_replace_all_regex
.add_products <- function(substrates, reaction = "RHEA:15421") {

    ## make the substrates argument a little bit smaller
    .s <- substrates
    ## start of if statements 
    ## depending on the reaction type, replace the expressions of the substrates
    ## and products in the reaction
    
    ## acyldhap_to_alkyldhap
    if (reaction %in% c("RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174")) { 
        .s$FA <- stringi::stri_replace_all_regex(str = .s$AcylDHAP, 
            pattern = "DHAP\\(", replacement = "FA(")
        .s$AlkylDHAP <- stringi::stri_replace_all_regex(str = .s$FAO, 
            pattern = "FAO\\(", replacement = "DHAP(O-")
    }

    ## alkyldhap_to_lpao
    if (reaction %in% c("RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178")) {
        .s$sn1LPAO <- stringi::stri_replace_all_fixed(str = .s$AlkylDHAP, 
            pattern = "DHAP", replacement = "PA")
        .s$sn1LPAO <- stringi::stri_replace_all_regex(str = .s$sn1LPAO, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }

    ## cerp_to_cer
    if (reaction %in% c("RHEA:33743", "RHEA:33744", "RHEA:33745", "RHEA:33746")) {
        .s$Cer <- stringi::stri_replace_all_fixed(str = .s$CerP, 
            pattern = "CerP", replacement = "Cer")
    }
    
    ## cdpdg_to_pgp
    if (reaction %in% c("RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596")) { 
        .s$PGP <- stringi::stri_replace_all_fixed(str = .s$CDPDG, 
            pattern = "CDP-DG", replacement = "PGP")
    }
    
    ## cdpdg_to_pi
    if (reaction %in% c("RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583")) { 
        .s$PI <- stringi::stri_replace_all_fixed(str = .s$CDPDG, 
            pattern = "CDP-DG", replacement = "PI")     
    }
    
    ## cer_to_cerp
    if (reaction %in% c("RHEA:17929", "RHEA:17930", "RHEA:17931", "RHEA:17932")) { 
        .s$CerP <- stringi::stri_replace_all_fixed(str = .s$Cer, 
            pattern = "Cer", replacement = "CerP")
    }
    
    ## cer_to_glccer
    if (reaction %in% c("RHEA:12088", "RHEA:12089", "RHEA:12090", "RHEA:12091")) {
        .s$GlcCer <- stringi::stri_replace_all_fixed(str = .s$Cer, 
            pattern = "Cer", replacement = "GlcCer")
    }
            
    ## cer_to_sm
    if (reaction %in% c("RHEA:18765", "RHEA:18766", "RHEA:18767", "RHEA:18768")) {
        .s$SM <- stringi::stri_replace_all_fixed(str = .s$Cer, 
            pattern = "Cer", replacement = "SM")
    }
    
    ## cl_to_lcl
    if (reaction %in% c("RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938")) {
        .s$LCL <- unlist(lapply(isolate_radyls(.s$CL), 
            function(f) paste0("CL(1'-[", f[1], "/", f[2], "],3'-[0:0/", f[4], "])")))
        .s$FA <- unlist(lapply(isolate_radyls(.s$CL), 
            function(f) paste0("FA(", f[3], ")")))
    }
    
    ## coa_to_acyldhap
    if (reaction %in% c("RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660")) {
        .s$AcylDHAP <- stringi::stri_replace_all_fixed(str = .s$AcylCoA, 
            pattern = "CoA", replacement = "DHAP")
    }
    
    ## coa_to_FAO
    if (reaction %in% c("RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719")) {
        .s$FAO <- stringi::stri_replace_all_fixed(str = .s$AcylCoA, 
            pattern = "CoA", replacement = "FAO")
    }
    
    ## coa_to_lpa
    if (reaction %in% c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328")) {
        .s$sn1LPA <- stringi::stri_replace_all_fixed(str = .s$AcylCoA, pattern = "CoA", 
            replacement = "PA")
        .s$sn1LPA <- stringi::stri_replace_all_regex(str = .s$sn1LPA, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }
    
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715",
            "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666")) {
        .s$sn1MG <- unlist(lapply(isolate_radyls(.s$DG), 
            function(f) {paste0("MG(", f[1], "/0:0/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$DG), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    ## dg_to_sn2mg
    if (reaction %in% c("RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278")) {
        .s$sn2MG <- unlist(lapply(isolate_radyls(.s$DG), 
            function(f) {paste0("MG(0:0/", f[2], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$DG), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    ## dg_to_pa
    if (reaction %in% c("RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275")) {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PA")
        .s$PA <- stringi::stri_replace_all_regex(str = .s$PA, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }
    
    ## dg_to_pc
    if (reaction %in% c("RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942")) {  
        .s$PC <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PC")
        .s$PC <- stringi::stri_replace_all_regex(str = .s$PC, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    ## dg_to_pe
    if (reaction %in% c("RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946")) { 
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "DG",
            replacement = "PE")
        .s$PE <- stringi::stri_replace_all_regex(str = .s$PE, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }
    
    ## dg_to_tg
    if (reaction %in% c("RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871")) { 
        .s$TG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "/0:0", 
            replacement = paste0("/", 
                unlist(isolate_radyls(.s$AcylCoA))))
         .s$TG <- stringi::stri_replace_all_regex(str = .s$TG, pattern = "^DG", 
            replacement = "TG")
    }
    
    ## dgo_to_pco
    if (reaction %in% c("RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182")) { 
        .s$PCO <- stringi::stri_replace_all_fixed(str = .s$DGO, pattern = "DG", 
            replacement = "PC")
        .s$PCO <- stringi::stri_replace_all_regex(str = .s$PCO, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    ## dgo_to_peo
    if (reaction %in% c("RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190")) { 
        .s$PEO <- stringi::stri_replace_all_fixed(str = .s$DGO, pattern = "DG", 
            replacement = "PE")
        .s$PEO <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }
    
    ## dhcer_to_cer
    if (reaction %in% c("RHEA:46544", "RHEA:46545", "RHEA:46546", "RHEA:46547")) { ######## ????
        .s$Cer <- stringi::stri_replace_all_fixed(str = .s$DhCer, 
                            pattern = ":0", replacement = ":1")
            ## pattern = "Cer\\(d16:0\\(3OH,4OH\\)\\(15Me\\)\\/",
            ## replacement = "Cer(d16:1(4E)(3OH,4OH)(15Me)/")
    }

    ## dhcer_to_dhsm
    if (reaction %in% c("RHEA:44620", "RHEA:44621", "RHEA:44622", "RHEA:44623")) { 
        .s$DhSM <- stringi::stri_replace_all_fixed(str = .s$DhCer, 
            pattern = "Cer", replacement = "SM")
    }
    
    ## dhsm_to_dhcer
    if (reaction %in% c("RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303")) { 
        .s$DhCer <- stringi::stri_replace_all_fixed(str = .s$DhSM, 
            pattern = "SM", replacement = "Cer")
    }

    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424",
            "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886")) {
        .s$AcylCoA <- stringi::stri_replace_all_fixed(str = .s$FA, "FA", 
            replacement = "CoA")
    }
    
    ## lcl_to_cl
    if (reaction %in% c("RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842")) {
        
        ## isolate core from LCL (used to assemble LC)
        .s$LCLs2 <- isolate_radyls(.s$LCL)
        
        ## isolate core from AcylCoA (used to assemble LC)
        .s$AcylCoAs2 <- isolate_radyls(.s$AcylCoA)
        
        ## create LC
        .s$CL <- unlist(lapply(seq_along(.s$LCLs2), 
            function(f) paste0(
                "CL(1'-[", .s$LCLs2[[f]][1], "/", .s$LCLs2[[f]][2], 
                "],3'-[", .s$AcylCoAs2[[f]], "/", .s$LCLs2[[f]][4], "])")))
        
        ## remove the LCLs2 and AcylCoAs2 helper entries
        .s$LCLs2 <- NULL
        .s$AcylCoAs2 <- NULL
    }
    
    ## lnape_to_gpnae
    if (reaction %in% c("RHEA:45420", "RHEA:45421", "RHEA:45422", "RHEA:45423")) { 
        .s$GPNAE <- unlist(lapply(isolate_radyls(.s$LNAPE), 
            function(f) {paste0("GPNAE(", f[3], ")")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$LNAPE), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
    
    ## lpa_to_pa
    if (reaction %in% c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712")) {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$sn1LPA, 
            pattern = "/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("/", f)})))
    }
    
    ## lpao_to_pao
    if (reaction %in% c("RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238")) { 
       .s$PAO <- stringi::stri_replace_all_fixed(str = .s$sn1LPAO,
           pattern = "/0:0",
           replacement = unlist(lapply(isolate_radyls(.s$AcylCoA),
               function(f) {paste0("/", f)})))
    }

    ## sn1lpc_to_fa
    if (reaction %in% c("RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180")) {
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn1LPC), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    ## sn2lpc_to_fa
    if (reaction %in% c("RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699")) {
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn2LPC), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    ## sn1lpc_to_pc
    if (reaction %in% c("RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940")) { 
        .s$PC <- stringi::stri_replace_all_fixed(str = .s$sn1LPC, 
            pattern = "/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("/", f)})))
    }

    ## sn1lpe_to_fa
    if (reaction %in% c("RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970")) {
        ## sn1 loss
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn1LPE), 
            function(f) {paste0("FA(", f[1], ")")}))
    }

    ## sn2lpe_to_fa
    if (reaction == "sn2lpe_to_fa") {
        ## sn2 loss 
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn2LPE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## sn1lpe_to_pe
    if (reaction %in% c("RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998")) { 
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$sn1LPE, 
            pattern = "/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA),
                function(f) {paste0("/", f)})))
    }

    ## sn1lpi_to_pi
    if (reaction %in% c("RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198")) {
        .s$PI <- stringi::stri_replace_all_regex(str = .s$sn1LPI,
            pattern = "/0:0\\)$",
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA),
                function(f) {paste0("/", f, ")")})))
    }
    
    ## lpeo_to_peo
    if (reaction == "lpeo_to_peo") {
        .s$PEO <- stringi::stri_replace_all_fixed(str = .s$sn1LPEO, 
            pattern = "/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("/", f)})))
    }

    ## lpep_to_pep
    if (reaction %in% c("RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248")) {
        .s$PEP <- stringi::stri_replace_all_fixed(str = .s$sn1LPEP, 
            pattern = "/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("/", f)})))
    }

    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466",
            "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946")) {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$sn1MG, 
            pattern = "/0:0/0:0", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("/", f, "/0:0")})))
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "MG",
            replacement = "DG")
    }
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950",
            "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744")) {
        .s$DG <- stringi::stri_replace_all_regex(str = .s$sn2MG, 
            pattern = "\\(0:0/", 
            replacement = unlist(lapply(isolate_radyls(.s$AcylCoA), 
                function(f) {paste0("(", f, "/")})))
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$DG, pattern = "MG", 
            replacement = "DG")
    }
    
    ## sn1mg_to_fa
    if (reaction %in% c("RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022")) {
        ## "sn1 loss"
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn1MG), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
    
    ## sn2mg_to_fa
    if (reaction %in% c("RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874")) {
        ## "sn2 loss"
        .s$FA <- unlist(lapply(isolate_radyls(.s$sn2MG), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## sn1mg_to_lpa
    if (reaction %in% c("RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750")) {
        .s$sn1LPA <- stringi::stri_replace_all_fixed(str = .s$sn1MG, 
            pattern = "MG", replacement = "PA")
        .s$sn1LPA <- stringi::stri_replace_all_regex(
            str = .s$sn1LPA, "/0:0\\)$", 
            replacement = "\\)")
    }

    ## sn2mg_to_sn1mg
    if (reaction == "sn2mg_to_sn1mg") {
        ## "sn2 loss"
        .s$sn1MG <- unlist(lapply(isolate_radyls(.s$sn2MG), 
            function(f) {paste0("MG(", f[2], "/0:0/0:0)")}))
    }

    ## nae_to_fa
    if (reaction %in% c("RHEA:17505", "RHEA:17506", "RHEA:17507", "RHEA:17508",
            "RHEA:39995", "RHEA:39996", "RHEA:39997", "RHEA:39998")) {
        .s$FA <- stringi::stri_replace_all_fixed(str = .s$NAE, pattern = "NAE", 
            replacement = "FA")
    }

    ## nape_to_lnape
    if (reaction %in% c("RHEA:45460", "RHEA:45461", "RHEA:45462", "RHEA:45463")) {
        .s$LNAPE <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("NAPE(", f[1], "/0:0/", f[3], ")")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## nape_to_nae
    if (reaction %in% c("RHEA:33159", "RHEA:33160", "RHEA:33161", "RHEA:33162")) {
        .s$NAE <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("NAE(", f[3], ")")}))
        .s$PA <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("PA(", f[1], "/", f[2], ")")}))
    }

    ## nape_to_pnae
    if (reaction == "nape_to_pnae") {
        .s$PNAE <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("PNAE(", f[3], ")")}))
        .s$DG <- unlist(lapply(isolate_radyls(.s$NAPE), 
            function(f) {paste0("DG(", f[1], "/", f[2], "/0:0)")}))
    }
    
    ## napeo_to_nae
    if (reaction == "napeo_to_nae") {
        .s$NAE <- unlist(lapply(isolate_radyls(.s$NAPEO), 
            function(f) {paste0("NAE(", f[3], ")")}))
        .s$PAO <- unlist(lapply(isolate_radyls(.s$NAPEO), 
            function(f) {paste0("PA(", f[1], "/", f[2], ")")}))
    }
    
    ## pa_to_cdpdg
    if (reaction %in% c("RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232")) {
        .s$CDPDG <- stringi::stri_replace_all_fixed(str = .s$PA, pattern = "PA",
            replacement = "CDP-DG")
    }

    ## pa_to_dg
    if (reaction %in% c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432")) {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PA, pattern = "PA", 
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    ## pao_to_dgo
    if (reaction %in% c("RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242")) {
        .s$DGO <- stringi::stri_replace_all_fixed(str = .s$PAO, pattern = "PA", 
            replacement = "DG")
        .s$DGO <- stringi::stri_replace_all_regex(str = .s$DGO, 
            pattern = "\\)$", replacement = "/0:0\\)")
    }
    
    ## pc_to_dg
    if (reaction %in% c("RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607")) { 
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC", 
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }
    
    ## pc_to_sn1lpc
    if (reaction %in% c("RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804")) { 
        ## "sn2 loss"
        .s$sn1LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## pc_to_sn2lpc
    if (reaction %in% c("RHEA:18689", "RHEA:18690", "RHEA:18691","RHEA:18692")) {
        ## "sn2 loss"
        .s$sn2LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("FA(", f[1], ")")}))
    }
    
    ## pc_to_pa
    if (reaction %in% c("RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448")) { 
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC", 
            replacement = "PA")
    }
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) { 
        .s$PS <- stringi::stri_replace_all_fixed(str = .s$PC, pattern = "PC",
            replacement = "PS")
    }
    
    ## pco_to_lpco
    if (reaction %in% c("RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234")) {
        ## "sn2 loss"
        .s$sn1LPCO <- unlist(lapply(isolate_radyls(.s$PCO), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PCO), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    ## lpco_to_lpao
    if (reaction %in% c("RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930")) {
        .s$sn1LPAO <- stringi::stri_replace_all_fixed(str = .s$sn1LPCO,
            pattern = "PC(O-", replacement = "PA(O-")   
    }
    
    ## lpco_to_mgo
    if (reaction %in% c("RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086")) {
        .s$sn1MGO <-  unlist(lapply(isolate_radyls(.s$sn1LPCO),
            function(f) {paste0("MG(", f[1], "/0:0/0:0)")}))
    }
    
    ## lpco_to_pco
    if (reaction %in% c("RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995")) {
        .s$PCO <- stringi::stri_replace_all_fixed(str = .s$sn1LPCO, 
            pattern = "/0:0)", 
            replacement = paste0("/", isolate_radyls(.s$AcylCoA), ")"))
    }
    
    ## pe_to_dg
    if (reaction %in% c("RHEA:78951", "RHEA:78952", "RHEA:78953", "RHEA:78954")) {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE",
            replacement = "DG")
        .s$DG <- stringi::stri_replace_all_regex(str = .s$DG, pattern = "\\)$",
            replacement = "/0:0\\)")
    }
    
    ## pe_to_sn1lpe
    if (reaction %in% c("RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607")) {
        ## "sn2" loss
        .s$sn1LPE <- unlist(lapply(isolate_radyls(.s$PE), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    ## pe_to_sn2lpe
    if (reaction %in% c("RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411")) {
        ## "sn1" loss
        .s$sn2LPE <- unlist(lapply(isolate_radyls(.s$PE), 
            function(f) {paste0("PE(0:0/", f[1], ")")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PE), 
            function(f) {paste0("FA(", f[2], ")")}))
    }
    
    ## pe_to_nape_sn1
    if (reaction %in% c("RHEA:45188", "RHEA:45189", "RHEA:45190", "RHEA:45191")) {
        .s$sn2LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPE <- stringi::stri_replace_all_regex(str = .s$PE, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPE <- stringi::stri_replace_all_fixed(str = .s$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }

    ## pe_to_nape_sn2
    if (reaction %in% c("RHEA:45192", "RHEA:45193", "RHEA:45194", "RHEA:45195")) {
        .s$sn1LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPE <- stringi::stri_replace_all_regex(str = .s$PE, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPE <- stringi::stri_replace_all_fixed(str = .s$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }
    
    ## pe_to_pa
    if (reaction == "pe_to_pa") {
        .s$PA <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE", 
            replacement = "PA")
    }
    
    ## pe_to_ps
    if (reaction %in% c("RHEA:27606", "RHEA:27607", "RHEA:27608", "RHEA:27609")) {
        .s$PS <- stringi::stri_replace_all_fixed(str = .s$PE, pattern = "PE", 
            replacement = "PS")
    }
    
    ## peo_to_lpeo
    if (reaction == "peo_to_lpeo") {
        ## "sn2 loss"
        .s$sn1LPEO <- unlist(lapply(isolate_radyls(.s$PEO), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PEO), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## peo_to_napeo_sn1
    if (reaction == "peo_to_napeo_sn1") {
        .s$sn2LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPEO <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPEO <- stringi::stri_replace_all_fixed(str = .s$NAPEO, 
            pattern = "PE", replacement = "NAPE")
    }
    
    ## peo_to_napeo_sn2
    if (reaction == "peo_to_napeo_sn2") {
        .s$sn1LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPEO <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPEO <- stringi::stri_replace_all_fixed(str = .s$NAPEO, 
            pattern = "PE", replacement = "NAPE")
    }
    
    ## peo_to_pep
    if (reaction %in% c("RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959")) {
        .s$PEP <- stringi::stri_replace_all_regex(str = .s$PEO, 
            pattern = "PE\\(O-", replacement = "PE(P-")
    }
    
    ## pep_to_lpep
    if (reaction %in% c("RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198")) {
        ## "sn2 loss"
        .s$sn1LPEP <- unlist(lapply(isolate_radyls(.s$PEP), 
            function(f) {paste0("PE(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PEP), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## lpep_to_fal
    if (reaction %in% c("RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908")) {
        .s$FAL <- unlist(lapply(isolate_radyls(.s$sn1LPEP),
            function(f) paste0("FAL(", f[1], ")")))
        .s$FAL <- stringi::stri_replace_all_fixed(str = .s$FAL,
            pattern = "P-", replacement = "")
    }

    ## lpep_to_lpap
    if (reaction %in% c("RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206")) {
        .s$sn1LPAP <- stringi::stri_replace_all_fixed(str = .s$sn1LPEP,
            pattern = "PE(P-", replacement = "PA(P-")
    }

    ## lpep_to_mgp
    if (reaction %in% c("RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202")) {
        .s$sn1MGP <- unlist(lapply(isolate_radyls(.s$sn1LPEP),
            function(f) {paste0("MG(", f[1], "/0:0/0:0)")}))
    }

    ## pep_to_napep_sn1
    if (reaction %in% c("RHEA:63596", "RHEA:63597", "RHEA:63598", "RHEA:63599")) {
        .s$sn2LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(0:0/", f[2], ")")}))
        .s$NAPEP <- stringi::stri_replace_all_regex(str = .s$PEP, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[1], ")")})))
        .s$NAPEP <- stringi::stri_replace_all_fixed(str = .s$NAPEP,
            pattern = "PE", replacement = "NAPE")
    }

    ## pep_to_napep_sn2
    if (reaction == "pep_to_napep_sn2") {
        .s$sn1LPC <- unlist(lapply(isolate_radyls(.s$PC), 
            function(f) {paste0("PC(", f[1], "/0:0)")}))
        .s$NAPEP <- stringi::stri_replace_all_regex(str = .s$PEP, 
            pattern = "\\)$", 
            replacement = unlist(lapply(isolate_radyls(.s$PC), 
                function(f) {paste0("/", f[2], ")")})))
        .s$NAPEP <- stringi::stri_replace_all_fixed(str = .s$NAPEP, 
            pattern = "PE", replacement = "NAPE")
    }

    ## pg_to_cl
    if (reaction %in% c("RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934")) {
        ## isolate core from PGs
        .s$PGs2 <- stringi::stri_replace_all_regex(str = .s$PG, 
            pattern = "^PG\\(", replacement = "")
        .s$PGs2 <- stringi::stri_replace_all_regex(str = .s$PGs2, 
            pattern = "\\)$", replacement = "")

        ## isolate core from CDPDGs
        .s$CDPDGs2 <- stringi::stri_replace_all_regex(str = .s$CDPDG, 
            pattern = "^CDP-DG\\(", replacement = "")
        .s$CDPDGs2 <- stringi::stri_replace_all_regex(str = .s$CDPDGs2, 
            pattern = "\\)$", replacement = "")

        ## create CL
        .s$CL <- paste0("CL(1'-[", .s$PGs2, "],3'-[", .s$CDPDGs2, "])")
    }

    ## pgp_to_pg
    if (reaction %in% c("RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754")) {
        .s$PG <- stringi::stri_replace_all_fixed(str = .s$PGP, pattern = "PGP", 
            replacement = "PG")
    }

    ## pi_to_dg
    if (reaction %in% c("RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487")) {
        .s$DG <- stringi::stri_replace_all_fixed(str = .s$PI, pattern = "PI",
            replacement = "DG")
    }

    ## pi_to_sn1lpi
    if (reaction %in% c("RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004")) {
        .s$sn1LPI <- unlist(lapply(isolate_radyls(.s$PI), 
            function(f) {paste0("PI(", f[1], "/0:0)")}))
        .s$FA <- unlist(lapply(isolate_radyls(.s$PI), 
            function(f) {paste0("FA(", f[2], ")")}))
    }

    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831")) {
        .s$PE <- stringi::stri_replace_all_fixed(str = .s$PS, pattern = "PS",
            replacement = "PE")
    }

    ## sm_to_cer
    if (reaction %in% c("RHEA:45644", "RHEA:45645", "RHEA:45646", "RHEA:45647")) {
        .s$Cer <- stringi::stri_replace_all_fixed(str = .s$SM, pattern = "SM",
            replacement = "Cer")
    }

    ## sphinga_to_dhcer
    if (reaction %in% c("RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427")) { ########## ???????
        .s$SPHs2 <- stringi::stri_replace_all_regex(str = .s$SPH, pattern = "SPH\\(", 
            replacement = "") |>
            stringi::stri_replace_all_regex(pattern = "\\)$", replacement = "") |>
            stringi::stri_replace_all_fixed(pattern = ":0", replacement = ":1")
        .s$AcylCoAs2 <- stringi::stri_replace_all_regex(str = .s$AcylCoA, 
            pattern = "CoA\\(", replacement = "") |>
            stringi::stri_replace_all_regex(pattern = "\\)$", replacement = "")
        .s$DhCer <- paste0("Cer(", .s$SPHs2, "/", .s$AcylCoAs2, ")")
            ##stringi::stri_replace_all_regex(str = .s$AcylCoA, 
            ##pattern = "CoA\\(", replacement = "Cer(16:0(3OH,4OH,15Me)/")
        
        ## remove the SPHs2 and AcylCoAs2 helper entries
        .s$SPHs2 <- NULL
        .s$AcylCoAs2 <- NULL
    }

    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274",
            "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867")) {
        ## "sn1 loss"
        .s$sn1Loss_DG <- unlist(lapply(isolate_radyls(.s$TG), 
            function(f) {paste0("DG(", f[3], "/", f[2], "/0:0)")}))
        .s$sn1Loss_FA <- unlist(lapply(isolate_radyls(.s$TG), 
            function(f) {paste0("FA(", f[1], ")")}))

        ## "sn3 loss"
        .s$sn3Loss_DG <- unlist(lapply(isolate_radyls(.s$TG), 
            function(f) {paste0("DG(", f[1], "/", f[2], "/0:0)")}))
        .s$sn3Loss_FA <- unlist(lapply(isolate_radyls(.s$TG), 
            function(f) {paste0("FA(", f[3], ")")}))
    }
    
    ## end of if statements
    ## return the data frame with added product(s)
    .s
}


