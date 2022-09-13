#' @name .create_substrates_combinations
#' 
#' @title Create substrates combinations
#' 
#' @description
#' The function returns a \code{data.frame} with all combinations from
#' \code{substrates}.
#' 
#' @details
#' The string replacement depend on the \code{reaction} argument.
#' 
#' The function will internally \code{expand.grid} to create the combinations.
#' 
#' @param substrates list of character vector(s)
#' @param constraints character vector with same length as \code{length(substrates)}
#' @param negate logical vector with same length as \code{length(substrates)}
#' 
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de} and 
#'     Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @return data.frame
#'
#' @examples 
#'  FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' 
#' ## create data.frame of substrates
#' LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
.create_substrates_combinations <- function(substrates = substrates, 
    constraints = character(length(substrates)), 
    negate = logical(length(substrates))) {
    
    ## check if constraints has correct lengths and adjust brackets format
    ## that it can be used by grep
    if (length(constraints) != length(substrates))
        stop("'length(constraints)' has to be the same as 'length(substrates)'")
    #constraints <- gsub(x = constraints, pattern = "[(]", replacement = "[(]")
    #constraints <- gsub(x = constraints, pattern = "[)]", replacement = "[)]")
    
    ## check if negate has correct length
    if (length(negate) != length(substrates)) 
        stop("'length(negate)' has to be the same as 'length(substrates)'")
    
    ## constrain variable substrates
    substrates_l <- lapply(seq_along(substrates), function(i)
        substrates[[i]][grep(x = substrates[[i]], 
            pattern = constraints[i], invert = negate[i])])
        
    ## perform combinatorics between the entries in substrates
    substrates_df <- expand.grid(substrates_l, stringsAsFactors = FALSE) |>
        data.frame(stringsAsFactors = FALSE)
    
    ## add colnames to substrates
    colnames(substrates_df) <- names(substrates)
    
    ## return the object
    substrates_df
}

#' @name .check_colnames_substrates_combinations
#' 
#' @title Check if correct colnames are in \code{substrates_combinations}
#' 
#' @description
#' The helper function checks if the correct columns are in \code{substrates} 
#' depending on the \code{reaction}.
#' 
#' @details 
#' The function will throw an error if \code{reaction} is not of length 1 and
#' is not an implemented method.
#' 
#' The function will invisibly return the correct columns for \code{reaction}.
#' 
#' \code{.check_colnames_substrates_combinations} is a helper function to 
#' test the integrity of \code{df}.
#' 
#' @param substrates \code{data.frame} or \code{list}
#' @param reaction \code{character(1)}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @return character containing the valid colnames
#'
#' @examples
#' FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
#'     "FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
#'     "FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")
#' substrates <- list(FA = FA)
#' 
#' ## create data.frame of substrates
#' df_substrates <- LipidNetworkPredictR:::.create_substrates_combinations(
#'     substrates = substrates, 
#'     constraints = "", negate = FALSE)
#' LipidNetworkPredictR:::.check_colnames_substrates_combinations(
#'     substrates = df_substrates, reaction = "RHEA:15421")
.check_colnames_substrates_combinations <- function(substrates, reaction = "RHEA:15421") {
    
    if (length(reaction) != 1)
        stop("'reaction' has to be of length 1")
    
    ## convert to data.frame if substrates is a list
    if (is.data.frame(substrates))
        substrates <- as.list(substrates)
        
    ## acyldhap_to_alkyldhap
    if (reaction %in% c("RHEA:36171", "RHEA:36172", "RHEA:36173", "RHEA:36174")) 
        cols <- c("AcylDHAP", "FAO")
    
    ## alkyldhap_to_lpao
    if (reaction %in% c("RHEA:36175", "RHEA:36176", "RHEA:36177", "RHEA:36178")) 
        cols <- c("AlkylDHAP")
    
    if (reaction == "c1p_to_cer")
        cols <- c("C1P")
    
    ## cdpdg_to_pgp
    if (reaction %in% c("RHEA:12593", "RHEA:12594", "RHEA:12595", "RHEA:12596")) 
        cols <- c("CDPDG")
    
    ## cdpdg_to_pi
    if (reaction %in% c("RHEA:11580", "RHEA:11581", "RHEA:11582", "RHEA:11583")) 
        cols <- c("CDPDG")
    
    if (reaction == "RHEA:17929") ## cer_to_c1p
        cols <- c("Cer")
    
    if (reaction == "RHEA:12088") ## cer_to_glccer
        cols <- c("Cer")
    
    if (reaction == "RHEA:18765") ## cer_to_sm
        cols <- c("Cer")
    
    ## cl_to_lcl
    if (reaction %in% c("RHEA:32935", "RHEA:32936", "RHEA:32937", "RHEA:32938"))
        cols <- c("CL")
    
    ## coa_to_acyldhap
    if (reaction %in% c("RHEA:17657", "RHEA:17658", "RHEA:17659", "RHEA:17660"))
        cols <- c("AcylCoA")
    
    ## coa_to_fao
    if (reaction %in% c("RHEA:52716", "RHEA:52717", "RHEA:52718", "RHEA:52719")) 
        cols <- c("AcylCoA")
    
    ## coa_to_lpa
    if (reaction %in% c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328")) 
        cols <- c("AcylCoA")
    
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715",
            "RHEA:35663", "RHEA:35664", "RHEA:35665", "RHEA:35666"))
        cols <- c("DG")
    
    ## dg_to_sn2mg
    if (reaction %in% c("RHEA:33275", "RHEA:33276", "RHEA:33277", "RHEA:33278"))
        cols <- c("DG")
    
    ## dg_to_pa
    if (reaction %in% c("RHEA:10272", "RHEA:10273", "RHEA:10274", "RHEA:10275"))
        cols <- c("DG")
    
    ## dg_to_pc
    if (reaction %in% c("RHEA:32939", "RHEA:32940", "RHEA:32941", "RHEA:32942")) 
        cols <- c("DG")
    
    ## dg_to_pe
    if (reaction %in% c("RHEA:32943", "RHEA:32944", "RHEA:32945", "RHEA:32946")) ## 
        cols <- c("DG")
    
    ## dg_to_tg
    if (reaction %in% c("RHEA:10868", "RHEA:10869", "RHEA:10870", "RHEA:10871")) 
        cols <- c("DG", "AcylCoA")
    
    ## dgo_to_pco
    if (reaction %in% c("RHEA:36179", "RHEA:36180", "RHEA:36181", "RHEA:36182"))
        cols <- c("DGO")
    
    ## dgo_to_peo
    if (reaction %in% c("RHEA:36187", "RHEA:36188", "RHEA:36189", "RHEA:36190"))
        cols <- c("DGO")
    
    if (reaction == "dhcer_to_cer")
        cols <- c("DhCer")
    
    if (reaction == "RHEA:44620") ## dhcer_to_dhsm
        cols <- c("DhCer")
    
    ## dhsm_to_dhcer
    if (reaction %in% c("RHEA:45300", "RHEA:45301", "RHEA:45302", "RHEA:45303"))
        cols <- c("DhSM")
    
    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424",
            "RHEA:38883", "RHEA:38884", "RHEA:38885", "RHEA:38886"))
        cols <- c("FA")
    
    ## lcl_to_cl
    if (reaction %in% c("RHEA:35839", "RHEA:35840", "RHEA:35841", "RHEA:35842")) {
        cols <- c("LCL", "AcylCoA")
    }
    
    if (reaction == "RHEA:45420") ## lnape_to_gpnae
        cols <- c("LNAPE")
    
    ## lpa_to_pa
    if (reaction %in% c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712")) 
        cols <- c("sn1LPA", "AcylCoA")
    
    ## lpao_to_pao
    if (reaction %in% c("RHEA:36235", "RHEA:36236", "RHEA:36237", "RHEA:36238")) 
        cols <- c("sn1LPAO", "AcylCoA")
    
    ## sn1lpc_to_fa
    if (reaction %in% c("RHEA:15177", "RHEA:15178", "RHEA:15179", "RHEA:15180"))
        cols <- c("sn1LPC")
    
    ## sn2lpc_to_fa
    if (reaction %in% c("RHEA:44696", "RHEA:44697", "RHEA:44698", "RHEA:44699"))
        cols <- c("sn2LPC")
    
    ## sn1lpc_to_pc
    if (reaction %in% c("RHEA:12937", "RHEA:12938", "RHEA:12939", "RHEA:12940")) 
        cols <- c("sn1LPC", "AcylCoA")
    
    ## sn1lpe_to_fa
    if (reaction %in% c("RHEA:32967", "RHEA:32968", "RHEA:32969", "RHEA:32970")) 
        cols <- c("sn1LPE")
    
    if (reaction == "sn2lpe_to_fa")
        cols <- c("sn2LPE")
    
    ## sn1lpe_to_pe
    if (reaction %in% c("RHEA:32995", "RHEA:32996", "RHEA:32997", "RHEA:32998"))
        cols <- c("sn1LPE", "AcylCoA")
    
    ## sn1lpi_to_pi
    if (reaction %in% c("RHEA:33195", "RHEA:33196", "RHEA:33197", "RHEA:33198"))
        cols <- c("sn1LPI", "AcylCoA")
    
    ## lpeo_to_peo
    if (reaction == "lpeo_to_peo")
        cols <- c("sn1LPEO", "AcylCoA")
    
    ## lpep_to_pep
    if (reaction %in% c("RHEA:16245", "RHEA:16246", "RHEA:16247", "RHEA:16248"))
        cols <- c("sn1LPEP", "AcylCoA")
    
    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466",
            "RHEA:39943", "RHEA:39944", "RHEA:39945", "RHEA:39946"))
        cols <- c("sn1MG", "AcylCoA")
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950",
            "RHEA:16741", "RHEA:16742", "RHEA:16743", "RHEA:16744"))
        cols <- c("sn2MG", "AcylCoA")
    
    ## sn1mg_to_fa
    if (reaction %in% c("RHEA:34019", "RHEA:34020", "RHEA:34021", "RHEA:34022"))
        cols <- c("sn1MG")
    
    ## sn2mg_to_fa
    if (reaction %in% c("RHEA:32871", "RHEA:32872", "RHEA:32873", "RHEA:32874"))
        cols <- c("sn2MG") 
    
    ## sn1mg_to_lpa
    if (reaction %in% c("RHEA:33747", "RHEA:33748", "RHEA:33749", "RHEA:33750"))
        cols <- c("sn1MG")
    
    if (reaction == "sn2mg_to_sn1mg")
        cols <- c("sn2MG")
    
    if (reaction == "nae_to_fa")
        cols <- c("NAE")
    
    if (reaction == "nape_to_lnape")
        cols <- c("NAPE")
    
    if (reaction == "nape_to_nae")
        cols <- c("NAPE")
    
    if (reaction == "nape_to_pnae")
        cols <- c("NAPE")
    
    if (reaction == "napeo_to_nae")
        cols <- c("NAPEO")
    
    ## pa_to_cdpdg
    if (reaction %in% c("RHEA:16229", "RHEA:16230", "RHEA:16231", "RHEA:16232"))
        cols <- c("PA")
    
    ## pa_to_dg
    if (reaction %in% c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"))
        cols <- c("PA")
    
    ## pao_to_dgo
    if (reaction %in% c("RHEA:36239", "RHEA:36240", "RHEA:36241", "RHEA:36242"))
        cols <- c("PAO")
    
    ## pc_to_dg
    if (reaction %in% c("RHEA:10604", "RHEA:10605", "RHEA:10606", "RHEA:10607"))
        cols <- c("PC")
    
    ## pc_to_sn1lpc
    if (reaction %in% c("RHEA:15801", "RHEA:15802", "RHEA:15803", "RHEA:15804")) 
        cols <- c("PC")
    
    ## pc_to_sn2lpc
    if (reaction %in% c("RHEA:18689", "RHEA:18690","RHEA:18691","RHEA:18692"))
        cols <- c("PC")
    
    ## pc_to_pa
    if (reaction %in% c("RHEA:14445", "RHEA:14446", "RHEA:14447", "RHEA:14448"))
        cols <- c("PC")
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) 
        cols <- c("PC")
    
    ## pco_to_lpco
    if (reaction %in% c("RHEA:36231", "RHEA:36232", "RHEA:36233", "RHEA:36234"))
        cols <- c("PCO")
    
    ## lpco_to_lpao
    if (reaction %in% c("RHEA:39927", "RHEA:39928", "RHEA:39929", "RHEA:39930"))
        cols <- c("sn1LPCO")
    
    ## lpco_to_mgo
    if (reaction %in% c("RHEA:36083", "RHEA:36084", "RHEA:36085", "RHEA:36086"))
        cols <- c("sn1LPCO")
    
    ## lpco_to_pco
    if (reaction %in% c("RHEA:23992", "RHEA:23993", "RHEA:23994", "RHEA:23995"))
        cols <- c("sn1LPCO", "AcylCoA")
    
    ## pe_to_dg
    if (reaction == "pe_to_dg")
        cols <- c("PE")
    
    ## pe_to_sn1lpe
    if (reaction %in% c("RHEA:44604", "RHEA:44605", "RHEA:44606", "RHEA:44607"))
        cols <- c("PE")
    
    ## pe_to_sn2lpe
    if (reaction %in% c("RHEA:44408", "RHEA:44409", "RHEA:44410", "RHEA:44411"))
        cols <- c("PE")
    
    ## pe_to_nape_sn1
    if (reaction == "pe_to_nape_sn1")
        cols <- c("PE", "PC")
    
    ## pe_to_nape_sn2
    if (reaction == "pe_to_nape_sn2")
        cols <- c("PE", "PC")
    
    ## pe_to_pa
    if (reaction == "pe_to_pa")
        cols <- c("PE")
    
    ## pe_to_ps
    if (reaction %in% c("RHEA:27606", "RHEA:27607", "RHEA:27608","RHEA:27609"))
        cols <- c("PE")
    
    ## peo_to_lpeo
    if (reaction == "peo_to_lpeo")
        cols <- c("PEO")
    
    ## peo_to_napeo_sn1
    if (reaction == "peo_to_napeo_sn1")
        cols <- c("PEO", "PC")
    
    ## peo_to_napeo_sn2
    if (reaction == "peo_to_napeo_sn2")
        cols <- c("PEO", "PC")
    
    ## peo_to_pep
    if (reaction %in% c("RHEA:22956", "RHEA:22957", "RHEA:22958", "RHEA:22959"))
        cols <- c("PEO")
    
    ## pep_to_lpep
    if (reaction %in% c("RHEA:36195", "RHEA:36196", "RHEA:36197", "RHEA:36198"))
        cols <- c("PEP")
    
    ## lpep_to_fal
    if (reaction %in% c("RHEA:16905", "RHEA:16906", "RHEA:16907", "RHEA:16908"))
        cols <- c("sn1LPEP")
    
    ## lpep_to_lpap
    if (reaction %in% c("RHEA:36203", "RHEA:36204", "RHEA:36205", "RHEA:36206"))
        cols <- c("sn1LPEP")
    
    ## lpep_to_mgp
    if (reaction %in% c("RHEA:36199", "RHEA:36200", "RHEA:36201", "RHEA:36202"))
        cols <- c("sn1LPEP")
    
    ## pep_to_napep_sn1
    if (reaction == "pep_to_napep_sn1")
        cols <- c("PEP", "PC")
    
    ## pep_to_napep_sn2
    if (reaction == "pep_to_napep_sn2")
        cols <- c("PEP", "PC")
    
    ## pg_to_cl
    if (reaction %in% c("RHEA:32931", "RHEA:32932", "RHEA:32933", "RHEA:32934"))
        cols <- c("PG", "CDPDG")
    
    ## pgp_to_pg
    if (reaction %in% c("RHEA:33751", "RHEA:33752", "RHEA:33753", "RHEA:33754"))
        cols <- c("PGP")
    
    ## pi_to_dg
    if (reaction %in% c("RHEA:43484", "RHEA:43485", "RHEA:43486", "RHEA:43487"))
        cols <- "PI"
    
    ## pi_to_sn1lpi
    if (reaction %in% c("RHEA:18001", "RHEA:18002", "RHEA:18003", "RHEA:18004"))
        cols <- c("PI")
    
    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"))
        cols <- c("PS")
    
    if (reaction == "sm_to_cer")
        cols <- c("SM")
    
    ## sphinga_to_dhcer
    if (reaction %in% c("RHEA:53424", "RHEA:53425", "RHEA:53426", "RHEA:53427"))
        cols <- c("AcylCoA")
    
    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274",
            "RHEA:44864", "RHEA:44865", "RHEA:44866", "RHEA:44867"))
        cols <- c("TG")
    
    if (!all(cols %in% names(substrates)))
        stop(paste("columns", paste(sprintf("%s", cols), collapse = ", "), 
            "not in 'colnames(substrates)'"))
    
    ## return cols
    invisible(cols)
}
