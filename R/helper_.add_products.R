
.add_products <- function(substrates, reaction = "fa_to_coa") {

    ## start of if statements 
    ## depending on the reaction treat the 
    if (reaction == "acdhap_to_alkyldhap") {
        substrates$FA <- gsub(x = substrates$ACDHAP, pattern = "DHAP\\(", 
            replacement = "FA(")
        substrates$ALKYLDHAP <- gsub(x = substrates$FATOH, pattern = "FATOH\\(", 
            replacement = "DHAP(O-")
    }

    if (reaction == "alkyldhap_to_lpao") {
        substrates$LPAO <- gsub(x = substrates$ALKYLDHAP, pattern = "DHAP", 
            replacement = "PA")
        substrates$LPAO <- gsub(x = substrates$LPAO, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "c1p_to_cer") {
        substrates$CER <- gsub(x = substrates$C1P, pattern = "C1P", 
            replacement = "Cer")
    }
    
    if (reaction == "cdpdg_to_pgp") {
        substrates$PGP <- gsub(substrates$CDPDG, "CDP-DG", replacement = "PGP")
    }
            
    if (reaction == "cdpdg_to_pi") {
        substrates$PI <- gsub(x = substrates$CDPDG, pattern = "CDP-DG",
            replacement = "PI")       
    }
            
    if (reaction == "cer_to_c1p") {
        substrates$C1P <- gsub(x = substrates$CER, "Cer", replacement = "C1P")
    }
            
    if (reaction == "cer_to_glccer") {
        substrates$GLCCER <- gsub(substrates$CER, "Cer", replacement = "GlcCer")
    }
            
    if (reaction == "cer_to_sm") {
        substrates$SM <- gsub(x = substrates$CER, "Cer", replacement = "SM")
    }
    
    if (reaction == "coa_to_acdhap") {
        substrates$DHAP <- gsub(x = substrates$CoA, pattern = "CoA",
            replacement = "DHAP")
    }
    
    if (reaction == "coa_to_fatoh") {
        substrates$FATOH <- gsub(x = substrates$CoA, pattern = "CoA", 
            replacement = "FATOH")
    }
    
    if (reaction == "coa_to_lpa") {
        substrates$LPA <- gsub(x = substrates$CoA, pattern = "CoA", 
            replacement = "PA")
        substrates$LPA <- gsub(x = substrates$LPA, "\\)$", 
            replacement = "/0:0\\)")
    }
            
    if (reaction == "dg_to_sn1mg") {
        substrates$MG <- unlist(lapply(substrates$DG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("MG(", fatty_acyl[1], "/0:0/0:0)")
        }))
        substrates$FA <- unlist(lapply(substrates$DG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")
        }))
    }
    
    if (reaction == "dg_to_sn2mg") {
        substrates$MG <- unlist(lapply(substrates$DG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("MG(0:0/", fatty_acyl[2], "/0:0)")
        }))
        substrates$FA <- unlist(lapply(substrates$DG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[1], ")")
        }))
    }
            
    if (reaction == "dg_to_pa") {
        substrates$PA <- gsub(x = substrates$DG, pattern = "DG",
            replacement = "PA")
        substrates$PA <- gsub(x = substrates$PA, "/0:0\\)$",
            replacement = "\\)")
    }

    if (reaction == "dg_to_pc") {
        substrates$PC <- gsub(x = substrates$DG, pattern = "DG",
            replacement = "PC")
        substrates$PC <- sub(substrates$PC,
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    if (reaction == "dg_to_pe") {
        substrates$PE <- gsub(x = substrates$DG, pattern = "DG",
            replacement = "PE")
        substrates$PE <- gsub(x = substrates$PE, pattern = "/0:0\\)$",
            replacement = "\\)")
    }
            
    if (reaction == "dg_to_tg") {
        substrates$TG <- gsub(x = substrates$DG, pattern = "/0:0", 
            replacement = paste0("/", unlist(lapply(substrates$CoA, 
                lipidomicsUtils::isolate_radyls))))
         substrates$TG <- gsub(x = substrates$TG, pattern = "^DG", 
            replacement = "TG")
    }
            
    if (reaction == "dgo_to_pco") {
        substrates$PCO <- gsub(x = substrates$DGO, pattern = "DG", 
            replacement = "PC")
        substrates$PCO <- stringr::str_replace(substrates$PCO, 
            pattern = "/0:0\\)$", replacement = "\\)")
    }

    if (reaction == "dgo_to_peo") {
        substrates$PEO <- gsub(x = substrates$DGO, pattern = "DG", 
            replacement = "PE")
        substrates$PEO <- gsub(x = substrates$PEO, pattern = "/0:0\\)$",
            replacement = "\\)")
    }
            
    if (reaction == "dhcer_to_cer") {
        substrates$CER <- gsub(x = substrates$DHCER, 
            pattern = "Cer\\(d16:0\\(3OH,4OH\\)\\(15Me\\)\\/",
            replacement = "Cer(d16:1(4E)(3OH,4OH)(15Me)/")
    }
            
    if (reaction == "dhcer_to_dhsm") {
        substrates$DHSM <- gsub(x = substrates$DHCER, pattern = "Cer", 
            replacement = "SM")
    }
            
        
    if (reaction == "dhsm_to_dhcer") {
        substrates$DHCER <- gsub(x = substrates$DHSM, pattern = "SM", 
            replacement = "Cer")
    }
        
    if (reaction == "fa_to_coa") {
        substrates$CoA <- gsub(x = substrates$FA, "FA", replacement = "CoA")
    }
    
    if (reaction == "lnape_to_gpnae") {
        substrates$GPNAE <- unlist(lapply(substrates$LNAPE, function(f) {
            paste0("GPNAE(", lipidomicsUtils::isolate_radyls(f)[3], ")")}))
        substrates$FA <- unlist(lapply(substrates$LNAPE, function(f) {
            paste0("FA(", lipidomicsUtils::isolate_radyls(f)[1], ")")}))   
    }
        
    if (reaction == "lpa_to_pa") {
        substrates$PA <- gsub(x = substrates$LPA, pattern = "/0:0", 
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }
            
    if (reaction == "lpao_to_pao") {
        substrates$PAO <- gsub(x = substrates$LPAO, pattern = "/0:0", 
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }

    if (reaction == "sn1lpc_to_fa") {
        substrates$FA <- unlist(lapply(substrates$sn1LPC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[1], ")")
        }))
    }

    if (reaction == "sn21pc_to_fa") {
        substrates$FA <- unlist(lapply(substrates$sn2LPC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")
        }))
    }

    if (reaction == "sn1lpc_to_pc") {
        substrates$PC <- gsub(x = substrates$sn1LPC, pattern = "/0:0",
            pattern = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }

    if (reaction == "sn1lpe_to_fa") {
        ## sn2 loss
        substrates$FA <- unlist(lapply(substrates$sn1LPE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[1], ")")
        }))
    }

    if (reaction == "sn2lpe_to_fa") {
        ## sn2 loss 
        substrates$FA <- unlist(lapply(substrates$sn2LPE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")
        }))
    }

    if (reaction == "sn1lpe_to_pe") {
        substrates$PE <- gsub(x = substrates$sn1LPE, pattern = "/0:0",
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }

    if (reaction == "lpeo_to_peo") {
        substrates$PEO <- gsub(x = substrates$LPEO, pattern = "/0:0",
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }

    if (reaction == "lpep_to_pep") {
        substrates$PEP <- gsub(x = substrates$LPEP, pattern = "/0:0",
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f))})))
    }

    if (reaction == "sn1mg_to_dg") {
        substrates$DG <- gsub(x = substrates$sn1MG, pattern = "/0:0/0:0", 
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f), "/0:0")})))
        substrates$DG <- gsub(x = substrates$DG, pattern = "MG",
            replacement = "DG")
    }
    
    if (reaction == "sn2mg_to_dg") {
        substrates$DG <- gsub(x = substrates$sn2MG, pattern = "\\(0:0/", 
            replacement = unlist(lapply(substrates$CoA, function(f) {
                paste0("(", lipidomicsUtils::isolate_radyls(f), "/")})))
        substrates$DG <- gsub(x = substrates$DG, pattern = "MG", 
            replacement = "DG")
    }
    
    if (reaction == "sn1mg_to_fa") {
        ## "sn2 loss"
        substrates$FA <- unlist(lapply(substrates$sn1MG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[1], ")")
        }))
    }

    if (reaction == "sn2mg_to_fa") {
        ## "sn2 loss"
        substrates$FA <- unlist(lapply(substrates$sn2MG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")
        }))
    }

    if (reaction == "sn1mg_to_lpa") {
        substrates$LPA <- gsub(x = substrates$sn1MG, pattern = "MG", 
            replacement = "PA")
        substrates$LPA <- gsub(x = substrates$LPA, "/0:0\\)$", 
            replacement = "\\)")
    }

    if (reaction == "sn2mg_to_sn1mg") {
        ## "sn2 loss"
        substrates$sn1MG <- unlist(lapply(substrates$sn2MG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("MG(", fatty_acyl[2], "/0:0/0:0)")
        }))
    }

    if (reaction == "nae_to_fa") {
        substrates$FA <- gsub(x = substrates$NAE, pattern = "NAE", 
            replacement = "FA")
    }

    if (reaction == "nape_to_lnape") {
        substrates$LNAPE <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("NAPE(", lipidomicsUtils::isolate_radyls(f)[1], "/0:0/",
                lipidomicsUtils::isolate_radyls(f)[3], ")")}))
        substrates$FA <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("FA(", lipidomicsUtils::isolate_radyls(f)[2], ")")}))
    }

    if (reaction == "nape_to_nae") {
        substrates$NAE <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("NAE(", lipidomicsUtils::isolate_radyls(f)[3], ")")}))
        substrates$PA <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("PA(", lipidomicsUtils::isolate_radyls(f)[1], "/",
                lipidomicsUtils::isolate_radyls(f)[2], ")")}))
    }

    if (reaction == "nape_to_pnae") {
        substrates$PNAE <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("PNAE(", lipidomicsUtils::isolate_radyls(f)[3], ")")}))
        substrates$DG <- unlist(lapply(substrates$NAPE, function(f) {
            paste0("DG(", lipidomicsUtils::isolate_radyls(f)[1], "/",
                lipidomicsUtils::isolate_radyls(f)[2], "/0:0)")}))
    }
    
    if (reaction == "napeo_to_nae") {
        substrates$NAE <- unlist(lapply(substrates$NAPEO, function(f) {
            paste0("NAE(", lipidomicsUtils::isolate_radyls(f)[3], ")")}))
        substrates$PAO <- unlist(lapply(substrates$NAPEO, function(f) {
            paste0("PA(", lipidomicsUtils::isolate_radyls(f)[1], "/",
                lipidomicsUtils::isolate_radyls(f)[2], ")")}))
    }
    
    if (reaction == "pa_to_cdpdg") {
        substrates$CDPDG <- gsub(x = substrates$PA, pattern = "PA",
            replacement = "CDP-DG")
    }

    if (reaction == "pa_to_dg") {
        substrates$DG <- gsub(x = substrates$PA, pattern = "PA", 
            replacement = "DG")
        substrates$DG <- gsub(x = substrates$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "pao_to_dgo") {
        substrates$DGO <- gsub(x = substrates$PAO, pattern = "PA", 
            replacement = "DG")
        substrates$DGO <- gsub(x = substrates$DGO, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "pc_to_dg") {
        substrates$DG <- gsub(x = substrates$PC, pattern = "PC", 
            replacement = "DG")
        substrates$DG <- sub(substrates$DG, pattern = "\\)$", 
            replacement = "/0:0\\)")
    }

    if (reaction == "pc_to_sn1lpc") {
        ## "sn2 loss"
        substrates$sn1LPC <- unlist(lapply(substrates$PC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PC(", fatty_acyl[1], "/0:0)")}))
        substrates$FA <- unlist(lapply(substrates$PC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyls[2], ")")}))
    }

    if (reaction == "pc_to_sn2lpc") {
        ## "sn2 loss"
        substrates$sn2LPC <- unlist(lapply(substrates$PC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PC(0:0/", fatty_acyl[2], ")")}))
        substrates$FA <- unlist(lapply(substrates$PC, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[1], ")")}))
    }
    
    if (reaction == "pc_to_pa") {
        substrates$PA <- gsub(x = substrates$PC, pattern = "PC", 
            replacement = "PA")
    }
    
    if (reaction == "pc_to_ps") {
        substrates$PS <- gsub(x = substrates$PC, pattern = "PC",
            replacement = "PS")
    }
    
    if (reaction == "pco_to_lpco") {
        ## "sn2 loss"
        substrates$LPCO <- unlist(lapply(substrates$PCO, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PC(", fatty_acyl[1], "/0:0)")}))
        substrates$FA <- unlist(lapply(substrates$PCO, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")}))
    }
    
    if (reaction == "pe_to_dg") {
        substrates$DG <- gsub(x = substrates$PE, pattern = "PE",
            replacement = "DG")
        substrates$DG <- sub(substrates$DG, pattern = "\\)$",
            replacement = "/0:0\\)")
    }
    
    if (reaction == "pe_to_sn1lpe") {
        ## "sn2" loss
        substrates$sn1LPE <- unlist(lapply(substrates$PE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PE(", fatty_acyl[1], "/0:0)")}))
        substrates$FA <- unlist(lapply(substrates$PE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")}))
    }

    if (reaction == "pe_to_sn2lpe") {
        substrates$sn2LPE <- unlist(lapply(Lsubstrates$PE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PE(0:0/", fatty_acyl[1], ")")}))
        substrates$FA <- unlist(lapply(substrates$PE, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyls[2], ")")}))
    }
    
    if (reaction == "pe_to_nape_sn1") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(0:0/", lipidomicsUtils::isolate_radyls(f)[2], ")")}))
        substrates$NAPE <- gsub(x = substrates$PE, pattern = "\\)$",
            unlist(lapply(substrates$PC, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f)[1], ")")})))
        substrates$NAPE <- gsub(x = substrates$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }

    if (reaction == "pe_to_nape_sn2") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(", lipidomicsUtils::isolate_radyls(f)[1], "/0:0)")}))
        substrates$NAPE <- gsub(x = substrates$PE, pattern = "\\)$", 
            replacement = unlist(lapply(substrates$PC, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f)[2], ")")})))
        substrates$NAPE <- gsub(x = substrates$NAPE, pattern = "PE", 
            replacement = "NAPE")
    }
    
    if (reaction == "pe_to_pa") {
        substrates$PA <- gsub(x = substrates$PE, pattern = "PE", 
            replacement = "PA")
    }
    
    if (reaction == "pe_to_ps") {
        substrates$PS <- gsub(x = substrates$PE, pattern = "PE", 
            replacement = "PS")
    }
    
    if (reaction == "peo_to_lpeo") {
        ## "sn2 loss"
        substrates$LPEO <- unlist(lapply(substrates$PEO, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PE(", fatty_acyl[1], "/0:0)")}))
        substrates$FA <- unlist(lapply(substrates$PEO, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acylf[2], ")")}))
    }

    if (reaction == "peo_to_napeo_sn1") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(0:0/", lipidomicsUtils::isolate_radyls(f)[2], ")")}))
        substrates$NAPEO <- gsub(x = substrates$PEO, pattern = "\\)$",
            replacement = unlist(lapply(substrates$PC, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f)[1], ")")})))
        substrates$NAPEO <- gsub(x = substrates$NAPEO, pattern = "PE", 
            replacement = "NAPE")
    }
    
    if (reaction == "peo_to_napeo_sn2") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(", lipidomicsUtils::isolate_radyls(f)[1], "/0:0)")}))
        substrates$NAPEO <- gsub(x = substrates$PEO, pattern = "\\)$",
            replacement = unlist(lapply(substrates$PC, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f)[2], ")")})))
        substrates$NAPEO <- gsub(x = substrates$NAPEO, pattern = "PE",
            replacement = "NAPE")
        
    }
    
    if (reaction == "peo_to_pep") {
        substrates$PEP <- gsub(x = substrates$PEO, pattern = "PE\\(O-", 
            replacement = "PE(P-")
    }
    
    if (reaction == "pep_to_lpep") {
        ## "sn2 loss"
        substrates$LPEP <- unlist(lapply(substrates$PEP, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("PE(", fatty_acyl[1], "/0:0)")}))
        substrates$FA <- unlist(lapply(substrates$PEP, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("FA(", fatty_acyl[2], ")")}))
    }

    if (reaction == "pep_to_napep_sn1") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(0:0/", lipidomicsUtils::isolate_radyls(f)[2], ")")}))
        substrates$NAPEP <- gsub(x = substrates$PEP, pattern = "\\)$",
            replacement = unlist(lapply(substrates$PC, function(f) {
                paste0("/",lipidomicsUtils::isolate_radyls(f)[1], ")")})))
        substrates$NAPEP <- gsub(x = substrates$NAPEP, pattern = "PE",
            replacement = "NAPE")
    }

    if (reaction == "pep_to_napep_sn2") {
        substrates$LPC <- unlist(lapply(substrates$PC, function(f) {
            paste0("PC(", lipidomicsUtils::isolate_radyls(f)[1], "/0:0)")}))
        substrates$NAPEP <- gsub(x = substrates$PEP, pattern = "\\)$",
            replacement = unlist(lapply(substrates$PC, function(f) {
                paste0("/", lipidomicsUtils::isolate_radyls(f)[2], ")")})))
        substrates$NAPEP <- gsub(x = substrates$NAPEP, pattern = "PE", 
            replacement = "NAPE")
    }
    
    if (reaction == "pg_to_cl") {
        ## isolate core from PGs
        substrates$PGs2 <- gsub(x = substrates$PG, pattern = "^PG\\(",
            replacement = "")
        substrates$PGs2 <- gsub(x = substrates$PGs2, pattern = "\\)$",
            replacement = "")
         
        ## isolate core from CDPDGs
        substrates$CDPDGs2 <- gsub(x = substrates$CDPDG, pattern = "^CDP-DG\\(",
            replacement = "")
        substrates$CDPDGs2 <- gsub(x = substrates$CDPDGs2, pattern = "\\)$",
            replacement = "")
        
        ## create CL
        substrates$CL <- paste0("CL(1'-[", substrates$PGs2, "],3'-[",
            substrates$CDPDGs2, pattern = "])")
    }
    
    if (reaction == "pgp_to_pg") {
        substrates$PG <- gsub(x = substrates$PGP, pattern = "PGP", 
            replacement = "PG")
    }

    if (reaction == "ps_to_pe") {
        substrates$PE <- gsub(x = substrates$PS, pattern = "PS",
            replacement = "PE")
    }
    
    if (reaction == "sm_to_cer") {
        substrates$CER <- gsub(x = substrates$SM, pattern = "SM",
            replacement = "Cer")
    }

    if (reaction == "sphinga_to_dhcer") {
        substrates$DHCER <- gsub(x = substrates$CoA, pattern = "CoA\\(",
            replacement = "Cer(d16:0(3OH,4OH)(15Me)/")
    }

    if (reaction == "tg_to_dg") {
        ## "sn1 loss"
        substrates$sn1Loss_dg <- unlist(lapply(substrates$TG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_radyls(f)
            paste0("DG(", fatty_acyls[3], "/", fatty_acyl[2], "/0:0)")}))
        substrates$sn1Loss_fa <- unlist(lapply(substrates$TG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_fatty_acyls(f)
            paste0("FA(", fatty_acyl[1], ")")}))
        
        ## "sn3 loss"
        substrates$sn3Loss_dg <- unlist(lapply(substrates$TG, function(f) {
            fatty_acyl <- lipidomicsUtils::isolate_fatty_acyls(TG)
            paste0("DG(", fatty_acyl[1], "/", fatty_acyls[2], "/0:0)")}))
        substrates$sn3Loss_fa <- unlist(lapply(substrates$TG, function(f) {
            fatty_acyls <- lipidomicsUtils::isolate_fatty_acyls(f)
            paste0("FA(", fatty_acyls[3], ")")}))
    }
    
    ## end of if statements
    ## return the data frame with added product(s)
    substrates
}


