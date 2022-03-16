## function create_reaction

test_that("create_reaction works", {
    
    ## acdhap_to_alkyldhap
    acdhap <- "DHAP(18:0)"
    fatoh <- "FATOH(16:0)"
    reaction_l <- create_reaction(substrates = list(ACDHAP = acdhap, FATOH = fatoh), 
        template = NULL, reaction = "acdhap_to_alkyldhap")
    expect_equal(reaction_l[[1]]$ALKYLDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "DHAP(18:0) + FATOH(16:0) <=> H+ + FA(18:0) + DHAP(O-16:0)")
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    reaction_l <- create_reaction(substrates = list(ALKYLDHAP = alkyldhap), 
        template = NULL, reaction =  "alkyldhap_to_lpao")
    expect_equal(reaction_l[[1]]$LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H+ + NADPH + DHAP(O-18:0) <=> PA(O-18:0/0:0) + NADP+")
    
    ## c1p_to_cer
    c1p <- "C1P(d16:0(3OH,4OH)(15Me)/18:0)"
    reaction_l <- create_reaction(substrates = list(C1P = c1p), 
        template = NULL, reaction =  "c1p_to_cer")
    expect_equal(reaction_l[[1]]$CER, "Cer(d16:0(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + C1P(d16:0(3OH,4OH)(15Me)/18:0) <=> Pi + Cer(d16:0(3OH,4OH)(15Me)/18:0)")
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- create_reaction(substrates = list(CDPDG = cdpdg), 
            template = NULL, reaction =  "cdpdg_to_pgp")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PGP(12:0(11Me)/14:0)")
    
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- create_reaction(substrates = list(CDPDG = cdpdg), template = NULL, 
        reaction =  "cdpdg_to_pi")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "Inositol + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PI(12:0(11Me)/14:0)") 
    
    ## cer_to_c1p
    cer <- "Cer(d16:0(3OH,4OH)(15Me)/18:0)"
    reaction_l <- create_reaction(substrates = list(CER = cer), 
        template = NULL, reaction =  "cer_to_c1p")
    expect_equal(reaction_l[[1]]$C1P, "C1P(d16:0(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "ATP + Cer(d16:0(3OH,4OH)(15Me)/18:0) <=> H+ + ADP + C1P(d16:0(3OH,4OH)(15Me)/18:0)")
    
    ## cer_to_glccer
    cer <- "Cer(d16:0(3OH,4OH)(15Me)/18:0)"
    reaction_l <- create_reaction(substrates = list(CER = cer), 
        template = NULL, reaction =  "cer_to_glccer")
    expect_equal(reaction_l[[1]]$GLCCER, "GlcCer(d16:0(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "UDP-Glucose + Cer(d16:0(3OH,4OH)(15Me)/18:0) <=> H+ + UDP + GlcCer(d16:0(3OH,4OH)(15Me)/18:0)")
    
    ## cer_to_sm
    cer <- "Cer(d16:0(3OH,4OH)(15Me)/18:0)"
    reaction_l <- create_reaction(substrates = list(CER = cer), 
        template = NULL, reaction =  "cer_to_sm")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:0(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula,
        "PC + Cer(d16:0(3OH,4OH)(15Me)/18:0) <=> DG + SM(d16:0(3OH,4OH)(15Me)/18:0)")
    
    ## coa_to_acdhap
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(CoA = coa), 
        template = NULL, reaction =  "coa_to_acdhap")
    expect_equal(reaction_l[[1]]$DHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "Dihydroxyacetone-P + CoA(18:0) <=> CoA + DHAP(18:0)")
    
    ## coa_to_fatoh
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(CoA = coa), 
        template = NULL, reaction =  "coa_to_fatoh")
    expect_equal(reaction_l[[1]]$FATOH, "FATOH(18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(18:0) + 2 NADPH + 2 H+ <=> FATOH(18:0) + 2 NADP+ + CoA")
    
    ## coa_to_lpa
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(CoA = coa), 
        template = NULL, reaction =  "coa_to_lpa")
    expect_equal(reaction_l[[1]]$LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "Glycerol-3-P + CoA(18:0) <=> CoA + PA(18:0/0:0)")
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)" 
    reaction_l <- create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "dg_to_sn1mg")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "dg_to_sn2mg")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "dg_to_pa")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "ATP + DG(18:0/16:0/0:0) <=> H+ + ADP + PA(18:0/16:0)")
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "dg_to_pc")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CDP-Chol + DG(18:0/16:0/0:0) <=> H+ + CMP + PC(18:0/16:0)")
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "dg_to_pe")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CDP-Ethn + DG(18:0/16:0/0:0) <=> H+ + CMP + PE(18:0/16:0)")
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    coa <- "CoA(14:0)"
    reaction_l <- create_reaction(substrates = list(DG = dg, CoA = coa), 
        template = NULL, reaction =  "dg_to_tg")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) <=> CoA + TG(18:0/16:0/14:0)")
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "dgo_to_pco")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula,
        "CDP-Chol + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PC(O-18:0/16:0)")
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "dgo_to_peo")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CDP-Ethn + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PE(O-18:0/16:0)")
    
    ## dhcer_to_cer
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- create_reaction(substrates = list(DHCER = dhcer), 
        template = NULL, reaction =  "dhcer_to_cer")
    expect_equal(reaction_l[[1]]$CER, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula,
        "H+ + NADH + O2 + Cer(d16:0(3OH,4OH)(15Me)/12:0) <=> 2 H2O + NAD+ + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- create_reaction(substrates = list(DHCER = dhcer),
        template = NULL, reaction =  "dhcer_to_dhsm")
    expect_equal(reaction_l[[1]]$DHSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) <=> DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## dhsm_to_dhcer
    dhsm <- "SM(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- create_reaction(substrates = list(DHSM = dhsm), 
        template = NULL, reaction =  "dhsm_to_dhcer")
    expect_equal(reaction_l[[1]]$DHCER, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + SM(d16:0(3OH,4OH)(15Me)/12:0) <=> P-Choline + H+ + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    reaction_l <- create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "fa_to_coa")
    expect_equal(reaction_l[[1]]$CoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "ATP + CoA + FA(18:0) <=> PPi + AMP + CoA(18:0)")
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    reaction_l <- create_reaction(substrates = list(LNAPE = lnape), 
        template = NULL, reaction =  "lnape_to_gpnae")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, ######### check
        "NAPE(14:0/0:0/0:0) + H2O <=> GPNAE(0:0) + FA(14:0)") 
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    coa <- "CoA(14:0)"
    reaction_l <- create_reaction(substrates = list(LPA = lpa, CoA = coa), 
        template = NULL, reaction =  "lpa_to_pa")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(14:0) + PA(18:0/0:0) <=> CoA + PA(18:0/14:0)")
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    coa <- "CoA(14:0)"
    reaction_l <- create_reaction(substrates = list(LPAO = lpao,  CoA = coa), 
        template = NULL, reaction =  "lpao_to_pao")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PA(O-18:0/0:0) + CoA(14:0) <=> PA(O-18:0/14:0) + CoA")
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = NULL, reaction =  "sn1lpc_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(14:0/0:0) <=> Glycerophosphoethanolamine + H+ + FA(14:0)")
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    reaction_l <- create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = NULL, reaction =  "sn2lpc_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(0:0/14:0) <=> Glycerophosphoethanolamine + H+ + FA(14:0)")
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(sn1LPC = sn1lpc, CoA = coa), 
        template = NULL, reaction =  "sn1lpc_to_pc")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PC(14:0/0:0) + CoA(18:0) <=> PC(14:0/18:0) + CoA")
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn1LPE = pe), 
        template = NULL, reaction =  "sn1lpe_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(14:0/0:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    ## sn2lpe_to_fa
    sn2lpe <- "PE(0:0/14:0)"
    reaction_l <- create_reaction(substrates = list(sn2LPE = sn2lpe), 
        template = NULL, reaction =  "sn2lpe_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(0:0/14:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    ## sn1lpe_to_pe
    sn1lpe <- "PE(14:0/0:0)"
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(sn1LPE = sn1lpe, CoA = coa), 
        template = NULL, reaction =  "sn1lpe_to_pe")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(18:0) + PE(14:0/0:0) <=> CoA + PE(14:0/18:0)")
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(LPEO = lpeo, CoA = coa), 
        template = NULL, reaction =  "lpeo_to_peo")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(18:0) + PE(O-16:0/0:0) <=> CoA + PE(O-16:0/18:0)")
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    coa <- "CoA(18:0)"
    reaction_l <- create_reaction(substrates = list(LPEP = lpep, CoA = coa), 
        template = NULL, reaction =  "lpep_to_pep")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    coa <- "CoA(16:0)"
    reaction_l <- create_reaction(substrates = list(sn1MG = sn1mg, CoA = coa), 
        template = NULL, reaction =  "sn1mg_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> CoA + DG(14:0/16:0/0:0)")
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    coa <- "CoA(16:0)"
    reaction_l <- create_reaction(substrates = list(sn2MG = sn2mg, CoA = coa), 
        template = NULL, reaction =  "sn2mg_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <=> CoA + DG(16:0/14:0/0:0)")  
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "sn1mg_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + MG(14:0/0:0/0:0) <=> Glycerol + H+ + FA(14:0)") 
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "sn2mg_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + MG(0:0/14:0/0:0) <=> Glycerol + H+ + FA(14:0)")
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "sn1mg_to_lpa")
    expect_equal(reaction_l[[1]]$LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "ATP + MG(14:0/0:0/0:0) <=> H+ + ADP + PA(14:0/0:0)")
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "sn2mg_to_sn1mg")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "MG(0:0/14:0/0:0) <=> MG(14:0/0:0/0:0)")
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    reaction_l <- create_reaction(substrates = list(NAE = nae), 
        template = NULL, reaction =  "nae_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + NAE(18:0) <=> Ethanolamine + H+ + FA(18:0)")
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_lnape")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAPE(14:0/0:0/18:0) + FA(16:0)")
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_nae")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAE(18:0) + PA(14:0/16:0)")
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_pnae")
    expect_equal(reaction_l[[1]]$PNAE, "PNAE(18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> PNAE(18:0) + DG(14:0/16:0/0:0)")
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(NAPEO = napeo), 
        template = NULL, reaction =  "napeo_to_nae")
    expect_equal(reaction_l[[1]]$NAE, "NAE(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "NAPE(O-18:0/16:0/14:0) + H2O <=> NAE(14:0) + PA(O-18:0/16:0)") 
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "pa_to_cdpdg")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CTP + PA(14:0/16:0) <=> PPi + CDP-DG(14:0/16:0)")
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "pa_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PA(14:0/16:0) <=> Pi + DG(14:0/16:0/0:0)")
    
    ## pao_to_dgo
    pao <- "PA(O-16:0/14:0)" #################### check
    reaction_l <- create_reaction(substrates = list(PAO = pao), 
        template = NULL, reaction = "pao_to_dgo")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PA(O-16:0/14:0) <=> Pi + DG(O-16:0/14:0/0:0)") 
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    reaction_l <- create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "pc_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(20:0/18:0) <=> P-Choline + DG(20:0/18:0/0:0)")
    
    ## pc_to_sn1lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "pc_to_sn1lpc")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(16:0/14:0) <=> PC(16:0/0:0) + FA(14:0)") 
    
    ## pc_to_sn2lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "pc_to_sn2lpc")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(16:0/14:0) <=> PC(0:0/14:0) + FA(16:0)")
    
    ## pc_to_pa
    pc <- "PC(16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "pc_to_pa")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(16:0/14:0) <=> Choline + PA(16:0/14:0)")
    
    ## pc_to_ps
    pc <- "PC(16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "pc_to_ps")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "L-Serine + PC(16:0/14:0) <=> Choline + PS(16:0/14:0)")  
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PCO = pco), 
        template = NULL, reaction = "pco_to_lpco")
    expect_equal(reaction_l[[1]]$LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PC(O-16:0/14:0) <=> H+ + PC(O-16:0/0:0) + FA(14:0)")
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(14:0/16:0) <=> P-Ethanolamine + DG(14:0/16:0/0:0)")
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_sn1lpe")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(14:0/0:0)")  
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_sn2lpe")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(0:0/14:0)")
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe, PC = pc), 
        template = NULL, reaction = "pe_to_nape_sn1") ####################
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/18:0) + PC(0:0/20:0)")
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe, PC = pc), 
        template = NULL, reaction = "pe_to_nape_sn2")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/20:0) + PC(18:0/0:0)")
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_pa")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula,
        "H2O + PE(14:0/16:0) <=> Ethanolamine + H+ + PA(14:0/16:0)")
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    reaction_l <- create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction =  "pe_to_ps")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "L-Serine + PE(14:0/16:0) <=> Ethanolamine + PS(14:0/16:0)")
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "peo_to_lpeo")
    expect_equal(reaction_l[[1]]$LPEO, "PE(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(O-16:0/14:0) <=> H+ + PE(O-16:0/0:0) + FA(14:0)")
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = NULL, reaction =  "peo_to_napeo_sn1") ################
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = NULL, reaction =  "peo_to_napeo_sn2")
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")    
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "peo_to_pep")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(O-16:0/14:0) + NADPH + H+ + O2 <=> PE(P-16:0/14:0) + NADP+ + 2 H2O")
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PEP = pep), 
        template = NULL, reaction =  "pep_to_lpep")
    expect_equal(reaction_l[[1]]$LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PE(P-16:0/14:0) <=> H+ + PE(P-16:0/0:0) + FA(14:0)")  
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = NULL, reaction =  "pep_to_napep_sn1") ####
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/20:0) + PC(0:0/18:0)")
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = NULL, reaction =  "pep_to_napep_sn2")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")  
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    reaction_l <- create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = NULL, reaction =  "pg_to_cl")
    ## check
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) <=> H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PGP = pgp), 
        template = NULL, reaction = "pgp_to_pg")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + PGP(16:0/14:0) <=> Pi + PG(16:0/14:0)")
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    reaction_l <- create_reaction(substrates = list(PS = ps), 
        template = NULL, reaction = "ps_to_pe")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H+ + PS(14:0/14:0) <=> CO2 + PE(14:0/14:0)")
    
    ## sm_to_cer
    sm <- "SM(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- create_reaction(substrates = list(SM = sm), 
        template = NULL, reaction = "sm_to_cer")
    expect_equal(reaction_l[[1]]$CER, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "H2O + SM(d16:0(3OH,4OH)(15Me)/12:0) <=> P-Choline + H+ + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## sphinga_to_dhcer
    coa <- "CoA(12:0)"
    reaction_l <- create_reaction(substrates = list(CoA = coa), 
        template = NULL, reaction = "sphinga_to_dhcer")
    expect_equal(reaction_l[[1]]$DHCER, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <=> CoA + Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    reaction_l <- create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "tg_to_dg")
    expect_equal(reaction_l[[1]]$sn1Loss_dg, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_fa, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_dg, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_fa, "FA(14:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula[1], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$ReactionFormula[2], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(14:0) + DG(18:0/16:0/0:0)")
})
