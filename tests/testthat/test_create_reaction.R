## function .create_reaction
test_that(".create_reaction works", {
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fatoh <- "FATOH(16:0)"
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FATOH = fatoh), 
        template = NULL, reaction = "RHEA:36171")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DHAP(18:0) + FATOH(16:0) = H+ + FA(18:0) + DHAP(O-16:0)")
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FATOH = fatoh), 
        template = NULL, reaction = "RHEA:36172")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DHAP(18:0) + FATOH(16:0) => H+ + FA(18:0) + DHAP(O-16:0)")
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FATOH = fatoh), 
        template = NULL, reaction = "RHEA:36173")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DHAP(18:0) + FATOH(16:0) <= H+ + FA(18:0) + DHAP(O-16:0)")
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FATOH = fatoh), 
        template = NULL, reaction = "RHEA:36174")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DHAP(18:0) + FATOH(16:0) <=> H+ + FA(18:0) + DHAP(O-16:0)")
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap), 
        template = NULL, reaction =  "RHEA:36175")
    expect_equal(reaction_l[[1]]$LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + NADPH + DHAP(O-18:0) = PA(O-18:0/0:0) + NADP+")
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap), 
        template = NULL, reaction =  "RHEA:36176")
    expect_equal(reaction_l[[1]]$LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + NADPH + DHAP(O-18:0) => PA(O-18:0/0:0) + NADP+")
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap), 
        template = NULL, reaction =  "RHEA:36177")
    expect_equal(reaction_l[[1]]$LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + NADPH + DHAP(O-18:0) <= PA(O-18:0/0:0) + NADP+")
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap), 
        template = NULL, reaction =  "RHEA:36178")
    expect_equal(reaction_l[[1]]$LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + NADPH + DHAP(O-18:0) <=> PA(O-18:0/0:0) + NADP+")
    
    ## c1p_to_cer
    c1p <- "C1P(16:0(3OH,4OH,15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(C1P = c1p), 
        template = NULL, reaction =  "c1p_to_cer")
    expect_equal(reaction_l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + C1P(16:0(3OH,4OH,15Me)/18:0) <=> Pi + Cer(16:0(3OH,4OH,15Me)/18:0)")
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:12593")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) = H+ + CMP + PGP(12:0(11Me)/14:0)")
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:12594")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) => H+ + CMP + PGP(12:0(11Me)/14:0)")
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:12595")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) <= H+ + CMP + PGP(12:0(11Me)/14:0)")
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:12596")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PGP(12:0(11Me)/14:0)")
    
    ## cdpdg_to_pi
    cdgdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:11580")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "myo-Inositol + CDP-DG(12:0(11Me)/14:0) = H+ + CMP + PI(12:0(11Me)/14:0)") 
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:11581")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "myo-Inositol + CDP-DG(12:0(11Me)/14:0) => H+ + CMP + PI(12:0(11Me)/14:0)") 
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:11582")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "myo-Inositol + CDP-DG(12:0(11Me)/14:0) <= H+ + CMP + PI(12:0(11Me)/14:0)") 
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg), 
        template = NULL, reaction =  "RHEA:11583")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "myo-Inositol + CDP-DG(12:0(11Me)/14:0) <=> H+ + CMP + PI(12:0(11Me)/14:0)") 
    
    ## cer_to_c1p
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer), 
        template = NULL, reaction =  "RHEA:17929")
    expect_equal(reaction_l[[1]]$C1P, "C1P(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + Cer(16:0(3OH,4OH,15Me)/18:0) <=> H+ + ADP + C1P(16:0(3OH,4OH,15Me)/18:0)")
    
    ## cer_to_glccer
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer), 
        template = NULL, reaction =  "RHEA:12088")
    expect_equal(reaction_l[[1]]$GLCCer, "GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "UDP-Glucose + Cer(16:0(3OH,4OH,15Me)/18:0) <=> H+ + UDP + GlcCer(16:0(3OH,4OH,15Me)/18:0)")
    
    ## cer_to_sm
    cer <- "Cer(16:0(3OH,4OH,15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer), 
        template = NULL, reaction =  "RHEA:18765")
    expect_equal(reaction_l[[1]]$SM, "SM(16:0(3OH,4OH,15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PC + Cer(16:0(3OH,4OH,15Me)/18:0) <=> DG + SM(16:0(3OH,4OH,15Me)/18:0)")
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    reaction_l <- .create_reaction(substrates = list(CL = cl), 
        template = NULL, reaction = "RHEA:32935")
    expect_equal(reaction_l[[1]]$LCL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    reaction_l <- .create_reaction(substrates = list(CL = cl), 
        template = NULL, reaction = "RHEA:32936")
    expect_equal(reaction_l[[1]]$LCL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    reaction_l <- .create_reaction(substrates = list(CL = cl), 
        template = NULL, reaction = "RHEA:32937")
    expect_equal(reaction_l[[1]]$LCL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    reaction_l <- .create_reaction(substrates = list(CL = cl), 
        template = NULL, reaction = "RHEA:32938")
    expect_equal(reaction_l[[1]]$LCL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + H+ + FA(18:4(6Z,9Z,12Z,15Z))")
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:17657")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) = CoA + DHAP(18:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:17658")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) => CoA + DHAP(18:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:17659")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) <= CoA + DHAP(18:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:17660")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Dihydroxyacetone-P + CoA(18:0) <=> CoA + DHAP(18:0)")
    
    ## coa_to_fatoh
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:52716")
    expect_equal(reaction_l[[1]]$FATOH, "FATOH(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ = FATOH(18:0) + 2 NADP+ + CoA")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:52717")
    expect_equal(reaction_l[[1]]$FATOH, "FATOH(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ => FATOH(18:0) + 2 NADP+ + CoA")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:52718")
    expect_equal(reaction_l[[1]]$FATOH, "FATOH(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CoA(18:0) + 2 NADPH + 2 H+ <= FATOH(18:0) + 2 NADP+ + CoA")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:52719")
    expect_equal(reaction_l[[1]]$FATOH, "FATOH(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 NADPH + 2 H+ <=> FATOH(18:0) + 2 NADP+ + CoA")
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:15325")
    expect_equal(reaction_l[[1]]$LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) = CoA + PA(18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:15326")
    expect_equal(reaction_l[[1]]$LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) => CoA + PA(18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:15327")
    expect_equal(reaction_l[[1]]$LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) <= CoA + PA(18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:15328")
    expect_equal(reaction_l[[1]]$LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Glycerol-3-P + CoA(18:0) <=> CoA + PA(18:0/0:0)")
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)" 
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:44712")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:35663")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:44713")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:35664")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:44714")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:35665")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:44715")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:35666")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(18:0/0:0/0:0) + FA(16:0)")
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:33275")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) = H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:33276")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) => H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:33277")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <= H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:33278")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + DG(18:0/16:0/0:0) <=> H+ + MG(0:0/16:0/0:0) + FA(18:0)")
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:10272")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) = H+ + ADP + PA(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:10273")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) => H+ + ADP + PA(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:10274")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) <= H+ + ADP + PA(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:10275")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + DG(18:0/16:0/0:0) <=> H+ + ADP + PA(18:0/16:0)")
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32939")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) = H+ + CMP + PC(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32940")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) => H+ + CMP + PC(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32941")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) <= H+ + CMP + PC(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32942")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Choline + DG(18:0/16:0/0:0) <=> H+ + CMP + PC(18:0/16:0)")
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32943")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) = H+ + CMP + PE(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32944")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) => H+ + CMP + PE(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32945")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) <= H+ + CMP + PE(18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = NULL, reaction =  "RHEA:32946")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(18:0/16:0/0:0) <=> H+ + CMP + PE(18:0/16:0)")
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:10868")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) = CoA + TG(18:0/16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:10869")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) => CoA + TG(18:0/16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:10870")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) <= CoA + TG(18:0/16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:10871")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(14:0) + DG(18:0/16:0/0:0) <=> CoA + TG(18:0/16:0/14:0)")
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "RHEA:36179")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CDP-Choline + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PC(O-18:0/16:0)")
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "RHEA:36187")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) = H+ + CMP + PE(O-18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "RHEA:36188")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) => H+ + CMP + PE(O-18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "RHEA:36189")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) <= H+ + CMP + PE(O-18:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = NULL, reaction =  "RHEA:36190")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-Ethanolamine + DG(O-18:0/16:0/0:0) <=> H+ + CMP + PE(O-18:0/16:0)")
    
    ## dhcer_to_cer 
    dhcer <- "Cer(16:0(3OH,4OH,15Me)/12:0)" ###################################################
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer), 
        template = NULL, reaction =  "dhcer_to_cer")
    expect_equal(reaction_l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H+ + NADH + O2 + Cer(16:0(3OH,4OH,15Me)/12:0) <=> 2 H2O + NAD+ + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(16:1(3OH,4OH,15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = NULL, reaction =  "RHEA:44620")
    expect_equal(reaction_l[[1]]$DhSM, "SM(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC + Cer(16:1(3OH,4OH,15Me)/12:0) <=> DG + SM(16:1(3OH,4OH,15Me)/12:0)")
    
    ## dhsm_to_dhcer
    dhsm <- "SM(16:1(3OH,4OH,15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = NULL, reaction =  "RHEA:45300")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) = P-Choline + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = NULL, reaction =  "RHEA:45301")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) => P-Choline + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = NULL, reaction =  "RHEA:45302")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) <= P-Choline + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = NULL, reaction =  "RHEA:45303")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:1(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(16:1(3OH,4OH,15Me)/12:0) <=> P-Choline + H+ + Cer(16:1(3OH,4OH,15Me)/12:0)")
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:15421")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) = PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:38883")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) = PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:15422")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) => PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:38884")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) => PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:15423")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) <= PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:38885")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) <= PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:15424")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) <=> PPi + AMP + CoA(18:0)")
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = NULL, reaction =  "RHEA:38886")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + CoA + FA(18:0) <=> PPi + AMP + CoA(18:0)")
    
    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])"
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    reaction_l <- .create_reaction(
        substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:35839")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA") 
    reaction_l <- .create_reaction(substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:35840")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA") 
    reaction_l <- .create_reaction(substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:35841")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA") 
    reaction_l <- .create_reaction(substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:35842")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA") 
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(LNAPE = lnape), 
        template = NULL, reaction =  "RHEA:45420")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "NAPE(14:0/0:0/0:0) + H2O <=> GPNAE(0:0) + FA(14:0)") 
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(substrates = list(LPA = lpa, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:19709")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) = CoA + PA(18:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(LPA = lpa, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:19710")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) => CoA + PA(18:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(LPA = lpa, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:19711")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <= CoA + PA(18:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(LPA = lpa, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:19712")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <=> CoA + PA(18:0/14:0)")
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(substrates = list(LPAO = lpao, AcylCoA = acylcoa),
       template = NULL, reaction =  "RHEA:36235")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
       "PA(O-18:0/0:0) + CoA(14:0) = PA(O-18:0/14:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(LPAO = lpao, AcylCoA = acylcoa),
        template = NULL, reaction =  "RHEA:36236")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PA(O-18:0/0:0) + CoA(14:0) => PA(O-18:0/14:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(LPAO = lpao,  AcylCoA = acylcoa),
        template = NULL, reaction =  "RHEA:36237")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PA(O-18:0/0:0) + CoA(14:0) <= PA(O-18:0/14:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(LPAO = lpao,  AcylCoA = acylcoa),
        template = NULL, reaction =  "RHEA:36238")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PA(O-18:0/0:0) + CoA(14:0) <=> PA(O-18:0/14:0) + CoA")
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = NULL, reaction =  "RHEA:15177")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(14:0/0:0) = Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = NULL, reaction =  "RHEA:15178")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(14:0/0:0) => Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = NULL, reaction =  "RHEA:15179")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(14:0/0:0) <= Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = NULL, reaction =  "RHEA:15180")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(14:0/0:0) <=> Glycerophosphocholine + H+ + FA(14:0)")
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = NULL, reaction =  "RHEA:44696")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(0:0/14:0) = Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = NULL, reaction =  "RHEA:44697")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(0:0/14:0) => Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = NULL, reaction =  "RHEA:44698")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(0:0/14:0) <= Glycerophosphocholine + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = NULL, reaction =  "RHEA:44699")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(0:0/14:0) <=> Glycerophosphocholine + H+ + FA(14:0)")
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:12937")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) = PC(14:0/18:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:12938")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) => PC(14:0/18:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:12939")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <= PC(14:0/18:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:12940")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <=> PC(14:0/18:0) + CoA")

    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = NULL, reaction =  "RHEA:32967")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/0:0) = H+ + FA(14:0) + Glycerophosphoethanolamine")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = NULL, reaction =  "RHEA:32968")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/0:0) => H+ + FA(14:0) + Glycerophosphoethanolamine")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = NULL, reaction =  "RHEA:32969")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/0:0) <= H+ + FA(14:0) + Glycerophosphoethanolamine")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe),
        template = NULL, reaction =  "RHEA:32970")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + PE(14:0/0:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    ## sn2lpe_to_fa
    sn2lpe <- "PE(0:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(sn2LPE = sn2lpe), 
        template = NULL, reaction =  "sn2lpe_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(0:0/14:0) <=> H+ + FA(14:0) + Glycerophosphoethanolamine")
    
    ## sn1lpe_to_pe
    sn1lpe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32995")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) = CoA + PE(14:0/18:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32996")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) => CoA + PE(14:0/18:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32997")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) <= CoA + PE(14:0/18:0)")
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32998")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(14:0/0:0) <=> CoA + PE(14:0/18:0)")
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:33195")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) = PI(16:0/18:1(9Z)) + CoA")
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:33196")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) => PI(16:0/18:1(9Z)) + CoA")
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:33197")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <= PI(16:0/18:1(9Z)) + CoA")
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:33198")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <=> PI(16:0/18:1(9Z)) + CoA")
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(
        substrates = list(LPEO = lpeo, AcylCoA = acylcoa), 
        template = NULL, reaction =  "lpeo_to_peo")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(O-16:0/0:0) <=> CoA + PE(O-16:0/18:0)")
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(LPEP = lpep, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16245")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    reaction_l <- .create_reaction(
        substrates = list(LPEP = lpep, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16246")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    reaction_l <- .create_reaction(
        substrates = list(LPEP = lpep, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16247")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    reaction_l <- .create_reaction(
        substrates = list(LPEP = lpep, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16248")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(P-16:0/0:0) <=> CoA + PE(P-16:0/18:0)")
    
    
    c("", "RHEA:16246", "RHEA:16247", "RHEA:16248")
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:38463")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:39943")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:38464")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:39944")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:38465")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:39945")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:38466")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> CoA + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:39946")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> CoA + DG(14:0/16:0/0:0)")
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32947")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) = CoA + DG(16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16741")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) = CoA + DG(16:0/14:0/0:0)")  
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32948")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) => CoA + DG(16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16742")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) => CoA + DG(16:0/14:0/0:0)")  
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32949")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <= CoA + DG(16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16743")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <= CoA + DG(16:0/14:0/0:0)")  
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:32950")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <=> CoA + DG(16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = NULL, reaction =  "RHEA:16744")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(16:0) + MG(0:0/14:0/0:0) <=> CoA + DG(16:0/14:0/0:0)")  
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:34019")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) = Glycerol + H+ + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:34020")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) => Glycerol + H+ + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:34021")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) <= Glycerol + H+ + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:34022")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(14:0/0:0/0:0) <=> Glycerol + H+ + FA(14:0)") 
    
    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "RHEA:32871")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) = Glycerol + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "RHEA:32872")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) => Glycerol + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "RHEA:32873")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <= Glycerol + H+ + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "RHEA:32874")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <=> Glycerol + H+ + FA(14:0)")
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:33747")
    expect_equal(reaction_l[[1]]$LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) = H+ + ADP + PA(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:33748")
    expect_equal(reaction_l[[1]]$LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) => H+ + ADP + PA(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:33749")
    expect_equal(reaction_l[[1]]$LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) <= H+ + ADP + PA(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = NULL, reaction =  "RHEA:33750")
    expect_equal(reaction_l[[1]]$LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "ATP + MG(14:0/0:0/0:0) <=> H+ + ADP + PA(14:0/0:0)")
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = NULL, reaction =  "sn2mg_to_sn1mg")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) <=> MG(14:0/0:0/0:0)")
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = NULL, reaction =  "nae_to_fa")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) <=> Ethanolamine + H+ + FA(18:0)")
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_lnape")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAPE(14:0/0:0/18:0) + FA(16:0)")
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_nae")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> NAE(18:0) + PA(14:0/16:0)")
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = NULL, reaction =  "nape_to_pnae")
    expect_equal(reaction_l[[1]]$PNAE, "PNAE(18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> PNAE(18:0) + DG(14:0/16:0/0:0)")
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(NAPEO = napeo), 
        template = NULL, reaction =  "napeo_to_nae")
    expect_equal(reaction_l[[1]]$NAE, "NAE(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(O-18:0/16:0/14:0) + H2O <=> NAE(14:0) + PA(O-18:0/16:0)") 
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:16229")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CTP + PA(14:0/16:0) = PPi + CDP-DG(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:16230")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CTP + PA(14:0/16:0) => PPi + CDP-DG(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:16231")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CTP + PA(14:0/16:0) <= PPi + CDP-DG(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:16232")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CTP + PA(14:0/16:0) <=> PPi + CDP-DG(14:0/16:0)")
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:27429")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(14:0/16:0) = Pi + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:27430")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(14:0/16:0) => Pi + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:27431")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(14:0/16:0) <= Pi + DG(14:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = NULL, reaction = "RHEA:27432")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(14:0/16:0) <=> Pi + DG(14:0/16:0/0:0)")
    
    ## pao_to_dgo
    pao <- "PA(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = NULL, reaction = "RHEA:36239")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(O-16:0/14:0) = Pi + DG(O-16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = NULL, reaction = "RHEA:36240")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(O-16:0/14:0) => Pi + DG(O-16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = NULL, reaction = "RHEA:36241")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(O-16:0/14:0) <= Pi + DG(O-16:0/14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = NULL, reaction = "RHEA:36242")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PA(O-16:0/14:0) <=> Pi + DG(O-16:0/14:0/0:0)") 
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:10604")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(20:0/18:0) = P-Choline + DG(20:0/18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:10605")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(20:0/18:0) => P-Choline + DG(20:0/18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:10606")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(20:0/18:0) <= P-Choline + DG(20:0/18:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:10607")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(20:0/18:0) <=> P-Choline + DG(20:0/18:0/0:0)")
    
    ## pc_to_sn1lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:15801")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) = PC(16:0/0:0) + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:15802")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) => PC(16:0/0:0) + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:15803")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <= PC(16:0/0:0) + FA(14:0)") 
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:15804")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <=> PC(16:0/0:0) + FA(14:0)") 
    
    ## pc_to_sn2lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:18689")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) = PC(0:0/14:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:18690")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) => PC(0:0/14:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:18691")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <= PC(0:0/14:0) + FA(16:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:18692")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <=> PC(0:0/14:0) + FA(16:0)")
    
    ## pc_to_pa
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:14445")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) = Choline + PA(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:14446")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) => Choline + PA(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:14447")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <= Choline + PA(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:14448")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(16:0/14:0) <=> Choline + PA(16:0/14:0)")
    
    ## pc_to_ps
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:45088")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PC(16:0/14:0) = Choline + PS(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:45089")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PC(16:0/14:0) => Choline + PS(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:45090")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PC(16:0/14:0) <= Choline + PS(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = NULL, reaction = "RHEA:45091")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PC(16:0/14:0) <=> Choline + PS(16:0/14:0)")
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = NULL, reaction = "RHEA:36231")
    expect_equal(reaction_l[[1]]$LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(O-16:0/14:0) = H+ + PC(O-16:0/0:0) + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = NULL, reaction = "RHEA:36232")
    expect_equal(reaction_l[[1]]$LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(O-16:0/14:0) => H+ + PC(O-16:0/0:0) + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = NULL, reaction = "RHEA:36233")
    expect_equal(reaction_l[[1]]$LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(O-16:0/14:0) <= H+ + PC(O-16:0/0:0) + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = NULL, reaction = "RHEA:36234")
    expect_equal(reaction_l[[1]]$LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PC(O-16:0/14:0) <=> H+ + PC(O-16:0/0:0) + FA(14:0)")
    
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_dg")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> P-Ethanolamine + DG(14:0/16:0/0:0)")
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44604")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) = H+ + FA(16:0) + PE(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44605")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) => H+ + FA(16:0) + PE(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44606")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <= H+ + FA(16:0) + PE(14:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44607")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(14:0/0:0)")
    
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44408")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) = H+ + FA(16:0) + PE(0:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44409")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) => H+ + FA(16:0) + PE(0:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44410")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <= H+ + FA(16:0) + PE(0:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:44411")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> H+ + FA(16:0) + PE(0:0/14:0)")
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = NULL, reaction = "pe_to_nape_sn1")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/18:0) + PC(0:0/20:0)")
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = NULL, reaction = "pe_to_nape_sn2")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + PC(18:0/20:0) <=> NAPE(14:0/16:0/20:0) + PC(18:0/0:0)")
    
    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "pe_to_pa")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + PE(14:0/16:0) <=> Ethanolamine + H+ + PA(14:0/16:0)")
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:27606")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PE(14:0/16:0) = Ethanolamine + PS(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:27607")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PE(14:0/16:0) => Ethanolamine + PS(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:27608")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PE(14:0/16:0) <= Ethanolamine + PS(14:0/16:0)")
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = NULL, reaction = "RHEA:27609")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "L-Serine + PE(14:0/16:0) <=> Ethanolamine + PS(14:0/16:0)")
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "peo_to_lpeo")
    expect_equal(reaction_l[[1]]$LPEO, "PE(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(O-16:0/14:0) <=> H+ + PE(O-16:0/0:0) + FA(14:0)")
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = NULL, reaction =  "peo_to_napeo_sn1") ################
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = NULL, reaction =  "peo_to_napeo_sn2")
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")    
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "RHEA:22956")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 = PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "RHEA:22957")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 => PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "RHEA:22958")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <= PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = NULL, reaction =  "RHEA:22959")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <=> PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = NULL, reaction =  "RHEA:36195")
    expect_equal(reaction_l[[1]]$LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(P-16:0/14:0) = H+ + PE(P-16:0/0:0) + FA(14:0)")
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = NULL, reaction =  "RHEA:36196")
    expect_equal(reaction_l[[1]]$LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(P-16:0/14:0) => H+ + PE(P-16:0/0:0) + FA(14:0)")  
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = NULL, reaction =  "RHEA:36197")
    expect_equal(reaction_l[[1]]$LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(P-16:0/14:0) <= H+ + PE(P-16:0/0:0) + FA(14:0)")  
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = NULL, reaction =  "RHEA:36198")
    expect_equal(reaction_l[[1]]$LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(P-16:0/14:0) <=> H+ + PE(P-16:0/0:0) + FA(14:0)")  
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = NULL, reaction =  "pep_to_napep_sn1") ####
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/20:0) + PC(0:0/18:0)")
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = NULL, reaction =  "pep_to_napep_sn2")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")  

    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = NULL, reaction = "RHEA:32931")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) = H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = NULL, reaction = "RHEA:32932")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) => H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = NULL, reaction = "RHEA:32933")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) <= H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = NULL, reaction = "RHEA:32934")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) + PG(18:4(6Z,9Z,12Z,15Z)/14:0) <=> H+ + CMP + CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = NULL, reaction = "RHEA:33751")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PGP(16:0/14:0) = Pi + PG(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = NULL, reaction = "RHEA:33752")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PGP(16:0/14:0) => Pi + PG(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = NULL, reaction = "RHEA:33753")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PGP(16:0/14:0) <= Pi + PG(16:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = NULL, reaction = "RHEA:33754")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PGP(16:0/14:0) <=> Pi + PG(16:0/14:0)")
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:43484")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) = myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:43485")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) => myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:43486")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <= myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:43487")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <=> myo-Inositol-1-P + H+ + DG(16:0/18:1(9Z))")
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:18001")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) = PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:18002")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) => PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:18003")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <= PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = NULL, reaction = "RHEA:18004")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PI(16:0/18:1(9Z)) <=> PI(16:0/0:0) + H+ + FA(18:1(9Z))")
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = NULL, reaction = "RHEA:20828")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + PS(14:0/14:0) = CO2 + PE(14:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = NULL, reaction = "RHEA:20829")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + PS(14:0/14:0) => CO2 + PE(14:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = NULL, reaction = "RHEA:20830")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + PS(14:0/14:0) <= CO2 + PE(14:0/14:0)")
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = NULL, reaction = "RHEA:20831")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H+ + PS(14:0/14:0) <=> CO2 + PE(14:0/14:0)")
    
    ## sm_to_cer
    sm <- "SM(16:0(3OH,4OH,15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(SM = sm), 
        template = NULL, reaction = "sm_to_cer")
    expect_equal(reaction_l[[1]]$Cer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + SM(16:0(3OH,4OH,15Me)/12:0) <=> P-Choline + H+ + Cer(16:0(3OH,4OH,15Me)/12:0)")
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:53424")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) = CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:53425")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) => CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:53426")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <= CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa), 
        template = NULL, reaction = "RHEA:53427")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(16:0(3OH,4OH,15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <=> CoA + Cer(16:0(3OH,4OH,15Me)/12:0)")
    ## should be: "CoA(12:0) + SPH(16:0(1OH,3OH,15Me)) <=> CoA + Cer(16:0(3OH,4OH,15Me)/12:0)"
    
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:33271")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) = H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) = H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:44864")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) = H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) = H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:33272")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) => H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) => H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:44865")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) => H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) => H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:33273")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) <= H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) <= H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:44866")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) <= H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) <= H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:33274")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(14:0) + DG(18:0/16:0/0:0)")
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = NULL, reaction = "RHEA:44867")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[2], 
        "H2O + TG(18:0/16:0/14:0) <=> H+ + FA(14:0) + DG(18:0/16:0/0:0)")
})


## function create_reactions
test_that("create_reactions works.", {
    FA <- c("FA(16:0)", "FA(12:0)", "FA(14:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_LPA + M_AcylCoA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE))
    
    reactions <- data.frame(order = reactions[, 1], RHEA = reactions[, 2],
        reactions = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
    
    ## run the function
    reactions_l <- create_reactions(substrates = list(FA = FA), 
        reactions = reactions)
    
    expect_equal(length(reactions_l), 4)
    
    ## first entry
    expect_equal(names(reactions_l[[1]][[1]]), "AcylCoA")
    expect_equal(reactions_l[[1]][[1]]$AcylCoA, 
        c("CoA(16:0)", "CoA(12:0)", "CoA(14:0)"))
    expect_equal(colnames(reactions_l[[1]][[2]]), 
        c("reaction_name", "reaction_formula", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_substrate", "reaction_product"))
    expect_equal(nrow(reactions_l[[1]][[2]]), 3)
    expect_equal(reactions_l[[1]][[2]]$reaction_formula[1:3],
        c("ATP + CoA + FA(16:0) = PPi + AMP + CoA(16:0)",
            "ATP + CoA + FA(12:0) = PPi + AMP + CoA(12:0)",
            "ATP + CoA + FA(14:0) = PPi + AMP + CoA(14:0)"))
    expect_equal(reactions_l[[1]][[2]]$reaction_isReversible[1:3], c("", "", ""))
    expect_equal(reactions_l[[1]][[2]]$reaction_geneAssociation[1:3], 
        c("", "", ""))
    expect_equal(reactions_l[[1]][[2]]$reaction_pathway[1:3], c("", "", ""))
    expect_equal(reactions_l[[1]][[2]]$reaction_substrate[1:3],
        c("ATP + CoA + FA(16:0)", "ATP + CoA + FA(12:0)", "ATP + CoA + FA(14:0)"))
    expect_equal(reactions_l[[1]][[2]]$reaction_product[1:3],
        c("PPi + AMP + CoA(16:0)", "PPi + AMP + CoA(12:0)", "PPi + AMP + CoA(14:0)"))
    
    ## second entry
    expect_equal(names(reactions_l[[2]][[1]]), "LPA")
    expect_equal(reactions_l[[2]][[1]]$LPA, 
        c("PA(16:0/0:0)", "PA(12:0/0:0)", "PA(14:0/0:0)"))
    expect_equal(colnames(reactions_l[[2]][[2]]), 
        c("reaction_name", "reaction_formula", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_substrate", "reaction_product"))
    expect_equal(nrow(reactions_l[[2]][[2]]), 3)
    expect_equal(reactions_l[[2]][[2]]$reaction_formula[1:3],
        c("Glycerol-3-P + CoA(16:0) = CoA + PA(16:0/0:0)",
            "Glycerol-3-P + CoA(12:0) = CoA + PA(12:0/0:0)",
            "Glycerol-3-P + CoA(14:0) = CoA + PA(14:0/0:0)"))
    expect_equal(reactions_l[[2]][[2]]$reaction_isReversible[1:3], 
        c("", "", ""))
    expect_equal(reactions_l[[2]][[2]]$reaction_geneAssociation[1:3], 
        c("", "", ""))
    expect_equal(reactions_l[[2]][[2]]$reaction_pathway[1:3], c("", "", ""))
    expect_equal(reactions_l[[2]][[2]]$reaction_substrate[1:3],
        c("Glycerol-3-P + CoA(16:0)", "Glycerol-3-P + CoA(12:0)",
            "Glycerol-3-P + CoA(14:0)"))
    expect_equal(reactions_l[[2]][[2]]$reaction_product[1:3],
        c("CoA + PA(16:0/0:0)", "CoA + PA(12:0/0:0)", "CoA + PA(14:0/0:0)"))
    
    ## third entry
    expect_equal(names(reactions_l[[3]][[1]]), "PA")
    expect_equal(reactions_l[[3]][[1]]$PA, 
        c("PA(16:0/16:0)", "PA(12:0/16:0)", "PA(14:0/16:0)", "PA(16:0/12:0)",
            "PA(12:0/12:0)", "PA(14:0/12:0)", "PA(16:0/14:0)", "PA(12:0/14:0)",
            "PA(14:0/14:0)"))
    expect_equal(colnames(reactions_l[[3]][[2]]), 
        c("reaction_name", "reaction_formula", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_substrate", "reaction_product"))
    expect_equal(nrow(reactions_l[[3]][[2]]), 9)
    expect_equal(reactions_l[[3]][[2]]$reaction_formula[1:3],
        c("PA(16:0/0:0) + CoA(16:0) = CoA + PA(16:0/16:0)",
            "PA(12:0/0:0) + CoA(16:0) = CoA + PA(12:0/16:0)",
            "PA(14:0/0:0) + CoA(16:0) = CoA + PA(14:0/16:0)"))
    expect_equal(reactions_l[[3]][[2]]$reaction_isReversible[1:3], c("", "", ""))
    expect_equal(reactions_l[[3]][[2]]$reaction_geneAssociation[1:3], 
                 c("", "", ""))
    expect_equal(reactions_l[[3]][[2]]$reaction_pathway[1:3], c("", "", ""))
    expect_equal(reactions_l[[3]][[2]]$reaction_substrate[1:3],
        c("PA(16:0/0:0) + CoA(16:0)", "PA(12:0/0:0) + CoA(16:0)",
            "PA(14:0/0:0) + CoA(16:0)"))
    expect_equal(reactions_l[[3]][[2]]$reaction_product[1:3],
        c("CoA + PA(16:0/16:0)", "CoA + PA(12:0/16:0)", "CoA + PA(14:0/16:0)"))
    
    ## fourth entry
    expect_equal(names(reactions_l[[4]][[1]]), "DG")
    expect_equal(reactions_l[[4]][[1]]$DG, 
        c("DG(16:0/16:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/16:0/0:0)",
            "DG(16:0/12:0/0:0)", "DG(12:0/12:0/0:0)", "DG(14:0/12:0/0:0)",
            "DG(16:0/14:0/0:0)", "DG(12:0/14:0/0:0)", "DG(14:0/14:0/0:0)"))
    expect_equal(colnames(reactions_l[[4]][[2]]), 
        c("reaction_name", "reaction_formula", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_substrate", "reaction_product"))
    expect_equal(nrow(reactions_l[[4]][[2]]), 9)
    expect_equal(reactions_l[[4]][[2]]$reaction_formula[1:3],
        c("H2O + PA(16:0/16:0) = Pi + DG(16:0/16:0/0:0)", 
            "H2O + PA(12:0/16:0) = Pi + DG(12:0/16:0/0:0)",
            "H2O + PA(14:0/16:0) = Pi + DG(14:0/16:0/0:0)"))
    expect_equal(reactions_l[[4]][[2]]$reaction_isReversible[1:3], c("", "", ""))
    expect_equal(reactions_l[[4]][[2]]$reaction_geneAssociation[1:3], 
                 c("", "", ""))
    expect_equal(reactions_l[[4]][[2]]$reaction_pathway[1:3], c("", "", ""))
    expect_equal(reactions_l[[4]][[2]]$reaction_substrate[1:3],
        c("H2O + PA(16:0/16:0)", "H2O + PA(12:0/16:0)", "H2O + PA(14:0/16:0)"))
    expect_equal(reactions_l[[4]][[2]]$reaction_product[1:3],
        c("Pi + DG(16:0/16:0/0:0)", "Pi + DG(12:0/16:0/0:0)",
            "Pi + DG(14:0/16:0/0:0)"))
})

## function create_reactions
test_that("create_reaction_adjacency_matrix works.", {
    FA <- c("FA(16:0)", "FA(12:0)", "FA(14:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_LPA + M_AcylCoA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE))
    
    reactions <- data.frame(order = reactions[, 1], RHEA = reactions[, 2],
                            reactions = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
    reactions_l <- create_reactions(substrates = list(FA = FA), 
                                    reactions = reactions)
    
    ## run the function
    adj <- create_reaction_adjacency_matrix(reaction_l = reactions_l)
    expect_equal(dim(adj), c(34, 34))
    expect_equal(rownames(adj), colnames(adj))
    expect_equal(rownames(adj), 
        c("AMP", "ATP", "CoA", "CoA(12:0)", "CoA(14:0)", "CoA(16:0)", 
            "DG(12:0/12:0/0:0)", "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)",
            "DG(14:0/12:0/0:0)", "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)",
            "DG(16:0/12:0/0:0)", "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)",
            "FA(12:0)", "FA(14:0)", "FA(16:0)", "Glycerol-3-P", "H2O",
            "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)", "PA(12:0/16:0)",
            "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", "PA(14:0/16:0)",
            "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)", "PA(16:0/16:0)",
            "PPi", "Pi"))
    expect_equal(sum(adj), 78)
    expect_equal(as.vector(adj[1:5, 1:5]), 
        c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 
            1, 1, 0, 0))
    expect_equal(adj["PA(14:0/0:0)", "PA(14:0/16:0)"], 1)
    expect_equal(adj["PA(14:0/0:0)", "Pi"], 0)
    expect_equal(adj["PA(14:0/0:0)", "PPi"], 0)
    expect_equal(adj["PA(14:0/0:0)", "PA(14:0/16:0)"], 1)
    expect_equal(as.vector(adj["ATP", c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)")]), 
        c(1, 1, 1))
    expect_equal(as.vector(adj[c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)"), "ATP"]), 
        c(0, 0, 0))
    
})
