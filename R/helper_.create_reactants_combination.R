#' @name .create_substrates_combinations
#' 
#' @title 
#' @param substrates list of character vector
#' @param constraints character vector with same length as \code{length(substrates)}
#' @param negate logical vector with same length as \code{length(substrates)}
.create_substrates_combinations <- function(substrates = substrates, 
    constraints = character(length(substrates)), 
    negate = logical(length(substrates))) {
    
    ## check if constraints has correct lengths and adjust brackets format
    ## that it can be used by grep
    if (length(constraints) != length(substrates))
        stop("'length(constraints)' has to be the same as 'length(substrates)'")
    constraints <- gsub(x = constraints, pattern = "[(]", replacement = "[(]")
    constraints <- gsub(x = constraints, pattern = "[)]", replacement = "[)]")
    
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

#
.check_colnames_substrates_combinations <- function(df, reaction = "fa_to_coa") {
    if (reaction == "acdhap_to_alkyldhap")
        .cols <- c("ACDHAP", "FATOH")
    
    if (reaction == "alkyldhap_to_lpao")
        .cols <- c("ALKYLDHAP")
    
    if (reaction == "c1p_to_cer")
        .cols <- c("C1P")
    
    if (reaction == "cdpdg_to_pgp")
        .cols <- c("CDPDG")
    
    if (reaction == "cdpdg_to_pi")
        .cols <- c("CDPDG")
    
    if (reaction == "cer_to_c1p")
        .cols <- c("CER")
    
    if (reaction == "cer_to_glccer")
        .cols <- c("CER")
    
    if (reaction == "cer_to_sm")
        .cols <- c("CER")
    
    if (reaction == "coa_to_acdhap")
        .cols <- c("CoA")
    
    if (reaction == "coa_to_fatoh")
        .cols <- c("CoA")
    
    if (reaction == "coa_to_lpa")
        .cols <- c("CoA")
    
    if (reaction == "dg_to_sn1mg")
        .cols <- c("DG")

    if (reaction == "dg_to_sn2mg")
        .cols <- c("DG")
    
    if (reaction == "dg_to_pa")
        .cols <- c("DG")
    
    if (reaction == "dg_to_pc")
        .cols <- c("DG")
    
    if (reaction == "dg_to_pe")
        .cols <- c("DG")
    
    if (reaction == "dg_to_tg")
        .cols <- c("DG", "CoA")
    
    if (reaction == "dgo_to_pco")
        .cols <- c("DGO")
    
    if (reaction == "dgo_to_peo")
        .cols <- c("DGO")
    
    if (reaction == "dhcer_to_cer")
        .cols <- c("DHCER")
    
    if (reaction == "dhcer_to_dhsm")
        .cols <- c("DHCER")
    
    if (reaction == "dhsm_to_dhcer")
        .cols <- c("DHSM")
    
    if (reaction == "fa_to_coa")
        .cols <- c("FA")
    
    if (reaction == "lnape_to_gpnae")
        .cols <- c("LNAPE")
    
    if (reaction == "lpa_to_pa")
        .cols <- c("LPA", "CoA")
    
    if (reaction == "lpao_to_pao")
        .cols <- c("LPAO", "CoA")
    
    if (reaction == "sn1lpc_to_fa")
        .cols <- c("sn1LPC")
    
    if (reaction == "sn21pc_to_fa")
        .cols <- c("sn2LPC")
    
    if (reaction == "sn1lpc_to_pc")
        .cols <- c("sn1LPC", "CoA")
    
    if (reaction == "sn1lpe_to_fa")
        .cols <- c("sn1LPE")
    
    if (reaction == "sn2lpe_to_fa")
        .cols <- c("sn2LPE")
    
    if (reaction == "sn1lpe_to_pe")
        .cols <- c("sn1LPE", "CoA")
    
    if (reaction == "lpeo_to_peo")
        .cols <- c("LPEO", "CoA")
    
    if (reaction == "lpep_to_pep")
        .cols <- c("LPEP", "CoA")
    
    if (reaction == "sn1mg_to_dg")
        .cols <- c("sn1MG", "CoA")
    
    if (reaction == "sn2mg_to_dg")
        .cols <- c("sn2MG", "CoA")
    
    if (reaction == "sn1mg_to_fa")
        .cols <- c("sn1MG")
    
    if (reaction == "sn2mg_to_fa")
        .cols <- c("sn2MG") 
    
    if (reaction == "sn1mg_to_lpa")
        .cols <- c("sn1MG")
    
    if (reaction == "sn2mg_to_sn1mg")
        .cols <- c("sn2MG")
    
    if (reaction == "nae_to_fa")
        .cols <- c("NAE")
    
    if (reaction == "nape_to_lnape")
        .cols <- c("NAPE")
    
    if (reaction == "nape_to_nae")
        .cols <- c("NAPE")
    
    if (reaction == "nape_to_pnae")
        .cols <- c("NAPE")
    
    if (reaction == "napeo_to_nae")
        .cols <- c("NAPEO")
    
    if (reaction == "pa_to_cdpdg")
        .cols <- c("PA")
    
    if (reaction == "pa_to_dg")
        .cols <- c("PA")
    
    if (reaction == "pao_to_dgo")
        .cols <- c("PAO")
    
    if (reaction == "pc_to_dg")
        .cols <- c("PC")
    
    if (reaction == "pc_to_sn1lpc")
        .cols <- c("PC")
    
    if (reaction == "pc_to_sn2lpc")
        .cols <- c("PC")
    
    if (reaction == "pc_to_pa")
        .cols <- c("PC")
    
    if (reaction == "pc_to_ps")
        .cols <- c("PC")
    
    if (reaction == "pco_to_lpco")
        .cols <- c("PCO")
    
    if (reaction == "pe_to_dg")
        .cols <- c("PE")
    
    if (reaction == "pe_to_sn1lpe")
        .cols <- c("PE", "")
    
    if (reaction == "pe_to_sn2lpe")
        .cols <- c("PE", "")
    
    if (reaction == "pe_to_nape_sn1")
        .cols <- c("PE", "PC")
    
    if (reaction == "pe_to_nape_sn2")
        .cols <- c("PE", "PC")
    
    if (reaction == "pe_to_pa")
        .cols <- c("PE")
    
    if (reaction == "pe_to_ps")
        .cols <- c("PE")
    
    if (reaction == "peo_to_lpeo")
        .cols <- c("PEO")
    
    if (reaction == "peo_to_napeo_sn1")
        .cols <- c("PEO", "PC")
    
    if (reaction == "peo_to_napeo_sn2")
        .cols <- c("PEO", "PC")
    
    if (reaction == "peo_to_pep")
        .cols <- c("PEO")
    
    if (reaction == "pep_to_lpep")
        .cols <- c("PEP", "")
    
    if (reaction == "pep_to_napep_sn1")
        .cols <- c("PEP", "PC")
    
    if (reaction == "pep_to_napep_sn2")
        .cols <- c("PEP", "PC")
    
    if (reaction == "pg_to_cl")
        .cols <- c("PG", "CDPDH")
    
    if (reaction == "pgp_to_pg")
        .cols <- c("PGP")
    
    if (reaction == "ps_to_pe")
        .cols <- c("PS")
    
    if (reaction == "sm_to_cer")
        .cols <- c("SM")
    
    if (reaction == "sphinga_to_dhcer")
        .cols <- c("CoA")
    
    if (reaction == "tg_to_dg")
        .cols <- c("TG")
    
    if (!all(.cols %in% colnames(df)))
        stop(paste("columns", paste(sprintf("%s", .cols), collapse = ", "), "not in 'colnames(df)'"))
}