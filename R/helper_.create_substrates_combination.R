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
        
    
    if (reaction == "RHEA:36171") ## acdhap_to_alkyldhap
        cols <- c("ACDHAP", "FATOH")
    
    if (reaction == "RHEA:36175") ## alkyldhap_to_lpao
        cols <- c("ALKYLDHAP")
    
    if (reaction == "c1p_to_cer")
        cols <- c("C1P")
    
    if (reaction == "RHEA:12593") ## cdpdg_to_pgp
        cols <- c("CDPDG")
    
    if (reaction == "RHEA:11580") ## cdpdg_to_pi
        cols <- c("CDPDG")
    
    if (reaction == "RHEA:17929") ## cer_to_c1p
        cols <- c("CER")
    
    if (reaction == "RHEA:12088") ## cer_to_glccer
        cols <- c("CER")
    
    if (reaction == "RHEA:18765") ## cer_to_sm
        cols <- c("CER")
    
    if (reaction == "coa_to_acdhap")
        cols <- c("CoA")
    
    if (reaction == "RHEA:52716") ## coa_to_fatoh
        cols <- c("CoA")
    
    ## coa_to_lpa
    if (reaction %in% c("RHEA:15325", "RHEA:15326", "RHEA:15327", "RHEA:15328")) 
        cols <- c("CoA")
    
    ## dg_to_sn1mg
    if (reaction %in% c("RHEA:44712", "RHEA:44713", "RHEA:44714", "RHEA:44715"))
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
        cols <- c("DG", "CoA")
    
    if (reaction == "RHEA:36179") ## dgo_to_pco
        cols <- c("DGO")
    
    if (reaction == "RHEA:36187") ## dgo_to_peo
        cols <- c("DGO")
    
    if (reaction == "dhcer_to_cer")
        cols <- c("DHCER")
    
    if (reaction == "RHEA:44620") ## dhcer_to_dhsm
        cols <- c("DHCER")
    
    if (reaction == "RHEA:19253") ## dhsm_to_dhcer
        cols <- c("DHSM")
    
    ## fa_to_coa
    if (reaction %in% c("RHEA:15421", "RHEA:15422", "RHEA:15423", "RHEA:15424"))
        cols <- c("FA")
    
    if (reaction == "RHEA:45420") ## lnape_to_gpnae
        cols <- c("LNAPE")
    
    ## lpa_to_pa
    if (reaction %in% c("RHEA:19709", "RHEA:19710", "RHEA:19711", "RHEA:19712")) 
        cols <- c("LPA", "CoA")
    
    if (reaction == "RHEA:36235") ## lpao_to_pao
        cols <- c("LPAO", "CoA")
    
    if (reaction == "RHEA:15177") ## sn1lpc_to_fa
        cols <- c("sn1LPC")
    
    if (reaction == "RHEA:44696") ## sn2lpc_to_fa
        cols <- c("sn2LPC")
    
    if (reaction == "RHEA:12937") ## sn1lpc_to_pc
        cols <- c("sn1LPC", "CoA")
    
    if (reaction == "RHEA:32967") ## sn1lpe_to_fa
        cols <- c("sn1LPE")
    
    if (reaction == "sn2lpe_to_fa")
        cols <- c("sn2LPE")
    
    if (reaction == "RHEA:32995") ## sn1lpe_to_pe
        cols <- c("sn1LPE", "CoA")
    
    if (reaction == "lpeo_to_peo")
        cols <- c("LPEO", "CoA")
    
    if (reaction == "lpep_to_pep")
        cols <- c("LPEP", "CoA")
    
    ## sn1mg_to_dg
    if (reaction %in% c("RHEA:38463", "RHEA:38464", "RHEA:38465", "RHEA:38466"))
        cols <- c("sn1MG", "CoA")
    
    ## sn2mg_to_dg
    if (reaction %in% c("RHEA:32947", "RHEA:32948", "RHEA:32949", "RHEA:32950"))
        cols <- c("sn2MG", "CoA")
    
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
    
    if (reaction == "pa_to_cdpdg")
        cols <- c("PA")
    
    ## pa_to_dg
    if (reaction %in% c("RHEA:27429", "RHEA:27430", "RHEA:27431", "RHEA:27432"))
        cols <- c("PA")
    
    if (reaction == "pao_to_dgo")
        cols <- c("PAO")
    
    if (reaction == "RHEA:10604") ## pc_to_dg
        cols <- c("PC")
    
    if (reaction == "RHEA:15801") ## pc_to_sn1lpc
        cols <- c("PC")
    
    if (reaction == "RHEA:18689") ## pc_to_sn2lpc
        cols <- c("PC")
    
    if (reaction == "RHEA:14445") ## pc_to_pa
        cols <- c("PC")
    
    ## pc_to_ps
    if (reaction %in% c("RHEA:45088", "RHEA:45089", "RHEA:45090", "RHEA:45091")) 
        cols <- c("PC")
    
    if (reaction == "pco_to_lpco")
        cols <- c("PCO")
    
    if (reaction == "pe_to_dg")
        cols <- c("PE")
    
    if (reaction == "pe_to_sn1lpe")
        cols <- c("PE")
    
    if (reaction == "RHEA:44408") ## pe_to_sn2lpe
        cols <- c("PE")
    
    if (reaction == "pe_to_nape_sn1")
        cols <- c("PE", "PC")
    
    if (reaction == "pe_to_nape_sn2")
        cols <- c("PE", "PC")
    
    if (reaction == "pe_to_pa")
        cols <- c("PE")
    
    ## pe_to_ps
    if (reaction == "RHEA:27606")
        cols <- c("PE")
    
    if (reaction == "peo_to_lpeo")
        cols <- c("PEO")
    
    if (reaction == "peo_to_napeo_sn1")
        cols <- c("PEO", "PC")
    
    if (reaction == "peo_to_napeo_sn2")
        cols <- c("PEO", "PC")
    
    if (reaction == "peo_to_pep")
        cols <- c("PEO")
    
    if (reaction == "pep_to_lpep")
        cols <- c("PEP")
    
    if (reaction == "pep_to_napep_sn1")
        cols <- c("PEP", "PC")
    
    if (reaction == "pep_to_napep_sn2")
        cols <- c("PEP", "PC")
    
    if (reaction == "pg_to_cl")
        cols <- c("PG", "CDPDG")
    
    if (reaction == "pgp_to_pg")
        cols <- c("PGP")
    
    ## ps_to_pe
    if (reaction %in% c("RHEA:20828", "RHEA:20829", "RHEA:20830", "RHEA:20831"))
        cols <- c("PS")
    
    if (reaction == "sm_to_cer")
        cols <- c("SM")
    
    if (reaction == "sphinga_to_dhcer")
        cols <- c("CoA")
    
    ## tg_to_dg
    if (reaction %in% c("RHEA:33271", "RHEA:33272", "RHEA:33273", "RHEA:33274"))
        cols <- c("TG")
    
    if (!all(cols %in% names(substrates)))
        stop(paste("columns", paste(sprintf("%s", cols), collapse = ", "), 
            "not in 'colnames(substrates)'"))
    
    ## return cols
    invisible(cols)
}
