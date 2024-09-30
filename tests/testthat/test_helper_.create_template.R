## function .create_template
test_that(".create_template works", {
    
    ## acyldhap_to_alkyldhap
    reaction <- "RHEA:36171"
    template <- .create_template(template = NA, reaction = reaction) 
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylDHAP + M_FAO = M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36172"
    template <- .create_template(template = NA, reaction = reaction) 
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylDHAP + M_FAO => M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36173"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylDHAP + M_FAO <= M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    reaction <- "RHEA:36174"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylDHAP + M_FAO <=> M_H+ + M_FA + M_AlkylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylDHAP", "M_FAO"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_AlkylDHAP"))
    
    ## alkyldhap_to_lpao
    reaction <- "RHEA:36175"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP = M_LPA-O + M_NADP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(template$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36176"
    template <- .create_template(template = NA, reaction = reaction) 
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP => M_LPA-O + M_NADP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(template$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36177"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP <= M_LPA-O + M_NADP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(template$reaction_product, c("M_LPA-O", "M_NADP"))
    
    reaction <- "RHEA:36178"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_NADPH + M_AlkylDHAP <=> M_LPA-O + M_NADP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_NADPH", "M_AlkylDHAP"))
    expect_equal(template$reaction_product, c("M_LPA-O", "M_NADP"))
    
    ## cerp_to_cer
    reaction <- "cerp_to_cer"
    template <- .create_template(template = NA, reaction = reaction) 
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_CerP <=> M_Pi + M_Cer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_CerP"))
    expect_equal(template$reaction_product, c("M_Pi", "M_Cer"))
    
    ## cdpdg_to_pgp
    reaction <- "RHEA:12593"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_CDP-DG = M_H+ + M_CMP + M_PGP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12594"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_CDP-DG => M_H+ + M_CMP + M_PGP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12595"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_CDP-DG <= M_H+ + M_CMP + M_PGP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    reaction <- "RHEA:12596"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_CDP-DG <=> M_H+ + M_CMP + M_PGP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PGP"))
    
    ## cdpdg_to_pi
    reaction <- "RHEA:11580"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_myo-Inositol + M_CDP-DG = M_H+ + M_CMP + M_PI")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11581"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_myo-Inositol + M_CDP-DG => M_H+ + M_CMP + M_PI")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11582"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_myo-Inositol + M_CDP-DG <= M_H+ + M_CMP + M_PI")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    reaction <- "RHEA:11583"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_myo-Inositol + M_CDP-DG <=> M_H+ + M_CMP + M_PI")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_myo-Inositol", "M_CDP-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PI"))
    
    ## cer_to_cerp
    reaction <- "RHEA:17929"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_Cer <=> M_H+ + M_ADP + M_CerP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_Cer"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_CerP"))
    
    ## cer_to_glccer
    reaction <- "RHEA:12088"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_UDP-Glucose + M_Cer <=> M_H+ + M_UDP + M_GlcCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_UDP-Glucose", "M_Cer"))
    expect_equal(template$reaction_product, c("M_H+", "M_UDP", "M_GlcCer"))
    
    ## cer_to_sm
    reaction <-  "RHEA:18765"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PC + M_Cer <=> M_1,2-DG + M_SM")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PC", "M_Cer"))
    expect_equal(template$reaction_product, c("M_1,2-DG", "M_SM"))
    
    ## cl_to_lcl
    reaction <- "RHEA:32935"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_CL = M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(template$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32936"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_CL => M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(template$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32937"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_CL <= M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(template$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32938"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_CL <=> M_1,2,4-LCL + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_CL" ))
    expect_equal(template$reaction_product, c("M_1,2,4-LCL", "M_H+", "M_FA"))
    
    ## coa_to_acyldhap
    reaction <- "RHEA:17657"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA = M_CoA + M_AcylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17658"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA => M_CoA + M_AcylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17659"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA <= M_CoA + M_AcylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    reaction <- "RHEA:17660"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Dihydroxyacetone-P + M_AcylCoA <=> M_CoA + M_AcylDHAP")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Dihydroxyacetone-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_AcylDHAP"))
    
    ## coa_to_FAO
    reaction <- "RHEA:52716"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ = M_FAO + 2 M_NADP + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(template$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    reaction <- "RHEA:52717"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ => M_FAO + 2 M_NADP + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(template$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    reaction <- "RHEA:52718"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ <= M_FAO + 2 M_NADP + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(template$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" )) 
    
    reaction <- "RHEA:52719"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + 2 M_NADPH + 2 M_H+ <=> M_FAO + 2 M_NADP + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "2 M_NADPH", "2 M_H+"))
    expect_equal(template$reaction_product, c("M_FAO", "2 M_NADP", "M_CoA" ))
    
    ## coa_to_lpa
    reaction <- "RHEA:15325"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15326"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_AcylCoA => M_CoA + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15327"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_AcylCoA <= M_CoA + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_LPA"))
    
    reaction <- "RHEA:15328"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_Glycerol-3-P + M_AcylCoA <=> M_CoA + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_Glycerol-3-P", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_LPA"))
    
    ## dg_to_sn1mg
    reaction <- "RHEA:44712"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:44713"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA")) 
    
    reaction <- "RHEA:44714"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:44715"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35663"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35664"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35665"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    reaction <- "RHEA:35666"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_1-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_1-MG", "M_FA"))
    
    ## dg_to_sn2mg
    reaction <- "RHEA:33275"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG = M_H+ + M_2-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33276"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG => M_H+ + M_2-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33277"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <= M_H+ + M_2-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    reaction <- "RHEA:33278"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1,2-DG <=> M_H+ + M_2-MG + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_2-MG", "M_FA"))
    
    ## dg_to_pa
    reaction <- "RHEA:10272"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1,2-DG = M_H+ + M_ADP + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10273"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1,2-DG => M_H+ + M_ADP + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10274"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1,2-DG <= M_H+ + M_ADP + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    reaction <- "RHEA:10275"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1,2-DG <=> M_H+ + M_ADP + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_PA"))
    
    ## dg_to_pc
    reaction <- "RHEA:32939"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_1,2-DG = M_H+ + M_CMP + M_PC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32940"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_1,2-DG => M_H+ + M_CMP + M_PC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32941"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_1,2-DG <= M_H+ + M_CMP + M_PC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    reaction <- "RHEA:32942"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_1,2-DG <=> M_H+ + M_CMP + M_PC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC"))
    
    ## dg_to_pe
    reaction <- "RHEA:32943"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG = M_H+ + M_CMP + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32944"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG => M_H+ + M_CMP + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32945"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG <= M_H+ + M_CMP + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    reaction <- "RHEA:32946"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_1,2-DG <=> M_H+ + M_CMP + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE"))
    
    ## dg_to_tg
    reaction <- "RHEA:10868"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1,2-DG = M_CoA + M_TG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10869"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1,2-DG => M_CoA + M_TG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10870"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1,2-DG <= M_CoA + M_TG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_TG"))
    
    reaction <- "RHEA:10871"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1,2-DG <=> M_CoA + M_TG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1,2-DG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_TG"))
    
    ## dgo_to_pco
    reaction <- "RHEA:36179"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_DG-O = M_H+ + M_CMP + M_PC-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36180"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_DG-O => M_H+ + M_CMP + M_PC-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36181"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_DG-O <= M_H+ + M_CMP + M_PC-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    reaction <- "RHEA:36182"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Choline + M_DG-O <=> M_H+ + M_CMP + M_PC-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Choline", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PC-O"))
    
    ## dgo_to_peo
    reaction <- "RHEA:36187"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_DG-O = M_H+ + M_CMP + M_PE-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36188"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_DG-O => M_H+ + M_CMP + M_PE-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36189"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_DG-O <= M_H+ + M_CMP + M_PE-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    reaction <- "RHEA:36190"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-Ethanolamine + M_DG-O <=> M_H+ + M_CMP + M_PE-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-Ethanolamine", "M_DG-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_PE-O"))
    
    ## dhcer_to_cer
    reaction <- "dhcer_to_cer"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_NADH + M_O2 + M_DhCer <=> 2 M_H2O + M_NAD + M_Cer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_NADH", "M_O2", "M_DhCer"))
    expect_equal(template$reaction_product, c("2 M_H2O", "M_NAD", "M_Cer"))
    
    ## dhcer_to_dhsm
    reaction <- "RHEA:44620"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PC + M_DhCer <=> M_1,2-DG + M_DhSM")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PC", "M_DhCer"))
    expect_equal(template$reaction_product, c("M_1,2-DG", "M_DhSM"))
    
    ## dhsm_to_dhcer
    reaction <- "RHEA:45300"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_DhSM = M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45301"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_DhSM => M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45302"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_DhSM <= M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))
    
    reaction <- "RHEA:45303"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_DhSM <=> M_Phosphocholine + M_H+ + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_DhSM"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_H+", "M_DhCer"))    
    
    ## fa_to_coa
    reaction <- "RHEA:15421"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA")) 
    
    reaction <- "RHEA:15422"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:15423"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:15424"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38883"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38884"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA => M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38885"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA <= M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    reaction <- "RHEA:38886"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_CoA + M_FA <=> M_PPi + M_AMP + M_AcylCoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_CoA", "M_FA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_AMP", "M_AcylCoA"))
    
    ## lcl_to_cl
    reaction <- "RHEA:35839"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA = M_CL + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CL", "M_CoA"))
    
    reaction <- "RHEA:35840"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA => M_CL + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CL", "M_CoA"))
    
    reaction <- "RHEA:35841"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA <= M_CL + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CL", "M_CoA"))
    
    reaction <- "RHEA:35842"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_1,2,4-LCL + M_AcylCoA <=> M_CL + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1,2,4-LCL", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CL", "M_CoA"))
    
    ## lnape_to_gpnae
    reaction <- "RHEA:45420"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LNAPE + M_H2O <=> M_GPNAE + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LNAPE", "M_H2O"))
    expect_equal(template$reaction_product, c("M_GPNAE", "M_FA"))
    
    ## lpa_to_pa
    reaction <- "RHEA:19709"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA + M_AcylCoA = M_CoA + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19710"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA + M_AcylCoA => M_CoA + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19711"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA + M_AcylCoA <= M_CoA + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PA"))
    
    reaction <- "RHEA:19712"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA + M_AcylCoA <=> M_CoA + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PA"))
    
    ## lpao_to_pao
    reaction <- "RHEA:36235"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA-O + M_AcylCoA = M_PA-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PA-O", "M_CoA"))
    
    reaction <- "RHEA:36236"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA-O + M_AcylCoA => M_PA-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PA-O", "M_CoA"))
    
    reaction <- "RHEA:36237"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA-O + M_AcylCoA <= M_PA-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PA-O", "M_CoA"))
    
    reaction <- "RHEA:36238"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,  "M_LPA-O + M_AcylCoA <=> M_PA-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_LPA-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PA-O", "M_CoA"))
    
    ## sn1lpc_to_fa
    reaction <- "RHEA:15177"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC = M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15178"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC => M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15179"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC <= M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:15180"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    ## sn2lpc_to_fa
    reaction <- "RHEA:44696"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-LPC = M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44697"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-LPC => M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44698"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-LPC <= M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    reaction <- "RHEA:44699"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-LPC <=> M_Glycerophosphocholine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-LPC"))
    expect_equal(template$reaction_product, c("M_Glycerophosphocholine", "M_H+", "M_FA"))
    
    ## sn1lpc_to_pc
    reaction <- "RHEA:12937"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC + M_AcylCoA = M_PC + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC", "M_CoA"))
    
    reaction <- "RHEA:12938"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC + M_AcylCoA => M_PC + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC", "M_CoA"))
    
    reaction <- "RHEA:12939"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC + M_AcylCoA <= M_PC + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC", "M_CoA"))
    
    reaction <- "RHEA:12940"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC + M_AcylCoA <=> M_PC + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC", "M_CoA"))
    
    ## sn1lpe_to_fa
    reaction <- "RHEA:32967"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE = M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32968"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE => M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32969"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE <= M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:32970"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    ## sn2lpe_to_fa
    reaction <- "sn2lpe_to_fa"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-LPE <=> M_H+ + M_FA + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-LPE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_Glycerophosphoethanolamine"))
    
    ## sn1lpe_to_pe
    reaction <- "RHEA:32995"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-LPE = M_CoA + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32996"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-LPE => M_CoA + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32997"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-LPE <= M_CoA + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE"))
    
    reaction <- "RHEA:32998"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-LPE <=> M_CoA + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-LPE"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE"))
    
    ## sn1lpi_to_pi
    reaction <- "RHEA:33195"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPI + M_AcylCoA = M_PI + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PI", "M_CoA"))
    
    reaction <- "RHEA:33196"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPI + M_AcylCoA => M_PI + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PI", "M_CoA"))
    
    reaction <- "RHEA:33197"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPI + M_AcylCoA <= M_PI + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PI", "M_CoA"))
    
    reaction <- "RHEA:33198"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPI + M_AcylCoA <=> M_PI + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPI", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PI", "M_CoA"))
    
    ## lpeo_to_peo
    reaction <- "lpeo_to_peo"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_ak2lgpe <=> M_CoA + M_PE-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_ak2lgpe"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE-O"))
    
    ## lpep_to_pep
    reaction <- "RHEA:16245"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_LPE-P = M_CoA + M_PE-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16246"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_LPE-P => M_CoA + M_PE-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16247"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_LPE-P <= M_CoA + M_PE-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE-P"))
    
    reaction <- "RHEA:16248"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_LPE-P <=> M_CoA + M_PE-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_LPE-P"))
    expect_equal(template$reaction_product, c("M_CoA", "M_PE-P"))
    
    ## sn1mg_to_dg
    reaction <- "RHEA:38463"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38464"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38465"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:38466"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39943"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA = M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39944"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA => M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39945"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA <= M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:39946"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-MG + M_AcylCoA <=> M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-MG", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    ## sn2mg_to_dg
    reaction <- "RHEA:32947"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32948"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32949"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:32950"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG <=> M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16741"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG = M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16742"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG => M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16743"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG <= M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    reaction <- "RHEA:16744"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_1-MG <=> M_CoA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_CoA", "M_1,2-DG"))
    
    ## sn1mg_to_fa
    reaction <- "RHEA:34019"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-MG = M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34020"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-MG => M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34021"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-MG <= M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:34022"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-MG <=> M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    ## sn2mg_to_fa
    reaction <- "RHEA:32871"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-MG = M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32872"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-MG => M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32873"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-MG <= M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    reaction <- "RHEA:32874"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_2-MG <=> M_Glycerol + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_2-MG"))
    expect_equal(template$reaction_product, c("M_Glycerol", "M_H+", "M_FA"))
    
    ## sn1mg_to_lpa
    reaction <- "RHEA:33747"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1-MG = M_H+ + M_ADP + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33748"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1-MG => M_H+ + M_ADP + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33749"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1-MG <= M_H+ + M_ADP + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    reaction <- "RHEA:33750"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_ATP + M_1-MG <=> M_H+ + M_ADP + M_LPA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_ATP", "M_1-MG"))
    expect_equal(template$reaction_product, c("M_H+", "M_ADP", "M_LPA"))
    
    ## sn2mg_to_sn1mg
    reaction <- "sn2mg_to_sn1mg"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_2-MG <=> M_1-MG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_2-MG"))
    expect_equal(template$reaction_product, c("M_1-MG"))
    
    ## nae_to_fa
    reaction <- "nae_to_fa"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_NAE <=> M_Ethanolamine + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_NAE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_H+", "M_FA"))
    
    ## nape_to_lnape
    reaction <- "nape_to_lnape"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_NAPE + M_H2O <=> M_LNAPE + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(template$reaction_product, c("M_LNAPE", "M_FA"))
    
    ## nape_to_nae
    reaction <- "nape_to_nae"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_NAPE + M_H2O <=> M_NAE + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(template$reaction_product, c("M_NAE", "M_PA"))
    
    ## nape_to_pnae
    reaction <- "nape_to_pnae"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_NAPE + M_H2O <=> M_PNAE + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_NAPE", "M_H2O"))
    expect_equal(template$reaction_product, c("M_PNAE", "M_1,2-DG"))
    
    ## napeo_to_nae
    reaction <- "napeo_to_nae"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_NAPEO + M_H2O <=> M_NAE + M_PA-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_NAPEO", "M_H2O"))
    expect_equal(template$reaction_product, c("M_NAE", "M_PA-O"))
    
    ## pa_to_cdpdg
    reaction <- "RHEA:16229"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CTP + M_PA = M_PPi + M_CDP-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16230"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CTP + M_PA => M_PPi + M_CDP-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16231"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CTP + M_PA <= M_PPi + M_CDP-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    reaction <- "RHEA:16232"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CTP + M_PA <=> M_PPi + M_CDP-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CTP", "M_PA"))
    expect_equal(template$reaction_product, c("M_PPi", "M_CDP-DG"))
    
    ## pa_to_dg
    reaction <- "RHEA:27429"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA = M_Pi + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(template$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27430"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA => M_Pi + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(template$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27431"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA <= M_Pi + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(template$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    reaction <- "RHEA:27432"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA <=> M_Pi + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA"))
    expect_equal(template$reaction_product, c("M_Pi", "M_1,2-DG"))
    
    ## pao_to_dgo
    reaction <- "RHEA:36239"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA-O = M_Pi + M_DG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(template$reaction_product, c("M_Pi", "M_DG-O"))
    
    reaction <- "RHEA:36240"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA-O => M_Pi + M_DG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(template$reaction_product, c("M_Pi", "M_DG-O"))
    
    reaction <- "RHEA:36241"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA-O <= M_Pi + M_DG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(template$reaction_product, c("M_Pi", "M_DG-O"))
    
    reaction <- "RHEA:36242"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PA-O <=> M_Pi + M_DG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PA-O"))
    expect_equal(template$reaction_product, c("M_Pi", "M_DG-O"))
    
    ## pc_to_dg
    reaction <- "RHEA:10604"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC = M_Phosphocholine + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10605"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC => M_Phosphocholine + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10606"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <= M_Phosphocholine + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    reaction <- "RHEA:10607"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <=> M_Phosphocholine + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_1,2-DG"))
    
    ## pc_to_sn1lpc
    reaction <- "RHEA:15801"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC = M_1-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15802"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC => M_1-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15803"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <= M_1-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_1-LPC", "M_FA"))
    
    reaction <- "RHEA:15804"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <=> M_1-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_1-LPC", "M_FA"))
    
    ## pc_to_sn2lpc
    reaction <- "RHEA:18689"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC = M_2-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18690"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC => M_2-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18691"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <= M_2-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_2-LPC", "M_FA"))
    
    reaction <- "RHEA:18692"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <=> M_2-LPC + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_2-LPC", "M_FA"))
    
    ## pc_to_pa
    reaction <- "RHEA:14445"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC = M_Choline + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14446"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC => M_Choline + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14447"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <= M_Choline + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PA"))
    
    reaction <- "RHEA:14448"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC <=> M_Choline + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PA"))
    
    ## pc_to_ps
    reaction <- "RHEA:45088"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PC = M_Choline + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45089"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PC => M_Choline + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45090"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PC <= M_Choline + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PS"))
    
    reaction <- "RHEA:45091"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PC <=> M_Choline + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PC"))
    expect_equal(template$reaction_product, c("M_Choline", "M_PS"))
    
    ## pco_to_lpco
    reaction <- "RHEA:36231"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC-O = M_H+ + M_LPC-O + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36232"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC-O => M_H+ + M_LPC-O + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36233"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC-O <= M_H+ + M_LPC-O + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    reaction <- "RHEA:36234"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PC-O <=> M_H+ + M_LPC-O + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPC-O", "M_FA"))
    
    ## lpco_to_lpao
    reaction <- "RHEA:39927"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O = M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39928"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O => M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39929"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O <= M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    reaction <- "RHEA:39930"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O <=> M_1-LPA-O + M_H+ + M_Choline")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_1-LPA-O", "M_H+", "M_Choline"))
    
    ## lpco_to_mgo
    reaction <- "RHEA:36083"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O = M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36084"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O => M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36085"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O <= M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    reaction <- "RHEA:36086"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPC-O <=> M_H+ + M_Phosphocholine + M_1-MG-O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPC-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_Phosphocholine", "M_1-MG-O"))
    
    ## lpco_to_pco
    reaction <- "RHEA:23992"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC-O + M_AcylCoA = M_PC-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC-O", "M_CoA"))
    
    reaction <- "RHEA:23993"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC-O + M_AcylCoA => M_PC-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC-O", "M_CoA"))
    
    reaction <- "RHEA:23994"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC-O + M_AcylCoA <= M_PC-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC-O", "M_CoA"))
    
    reaction <- "RHEA:23995"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_1-LPC-O + M_AcylCoA <=> M_PC-O + M_CoA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_1-LPC-O", "M_AcylCoA"))
    expect_equal(template$reaction_product, c("M_PC-O", "M_CoA"))
    
    ## pe_to_dg
    reaction <- "pe_to_dg"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <=> M_P-Ethanolamine + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_P-Ethanolamine", "M_1,2-DG"))
    
    ## pe_to_sn1lpe
    reaction <- "RHEA:44604"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE = M_H+ + M_FA + M_1-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44605"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE => M_H+ + M_FA + M_1-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44606"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <= M_H+ + M_FA + M_1-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    reaction <- "RHEA:44607"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <=> M_H+ + M_FA + M_1-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1-LPE"))
    
    ## pe_to_sn2lpe
    reaction <- "RHEA:44408"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE = M_H+ + M_FA + M_2-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44409"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE => M_H+ + M_FA + M_2-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44410"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <= M_H+ + M_FA + M_2-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    reaction <- "RHEA:44411"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <=> M_H+ + M_FA + M_2-LPE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_2-LPE"))
    
    ## pe_to_nape_sn1
    reaction <- "pe_to_nape_sn1"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE + M_PC <=> M_NAPE + M_2-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPE", "M_2-LPC"))
    
    ## pe_to_nape_sn2
    reaction <- "pe_to_nape_sn2"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE + M_PC <=> M_NAPE + M_1-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPE", "M_1-LPC"))
    
    ## pe_to_pa
    reaction <- "pe_to_pa"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE <=> M_Ethanolamine + M_H+ + M_PA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_H+", "M_PA"))
    
    ## pe_to_ps
    reaction <- "RHEA:27606"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PE = M_Ethanolamine + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27607"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PE => M_Ethanolamine + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27608"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PE <= M_Ethanolamine + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    reaction <- "RHEA:27609"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_L-Serine + M_PE <=> M_Ethanolamine + M_PS")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_L-Serine", "M_PE"))
    expect_equal(template$reaction_product, c("M_Ethanolamine", "M_PS"))
    
    ## peo_to_lpeo
    reaction <- "peo_to_lpeo"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE-O <=> M_H+ + M_LPE-O + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE-O"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPE-O", "M_FA"))
    
    ## peo_to_napeo_sn1
    reaction <- "peo_to_napeo_sn1"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_PC <=> M_NAPEO + M_2-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPEO", "M_2-LPC"))
    
    ## peo_to_napeo_sn2
    reaction <- "peo_to_napeo_sn2"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_PC <=> M_NAPEO + M_1-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPEO", "M_1-LPC"))
    
    ## peo_to_pep
    reaction <- "RHEA:22956"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 = M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(template$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    
    reaction <- "RHEA:22957"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 => M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(template$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    
    reaction <- "RHEA:22958"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <= M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(template$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    
    reaction <- "RHEA:22959"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-O + M_Fe2+-cytochrome_b5 + 2 M_H+ + M_O2 <=> M_PE-P + M_Fe3+-cytochrome_b5 + 2 M_H2O")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-O", "M_Fe2+-cytochrome_b5", "2 M_H+", "M_O2"))
    expect_equal(template$reaction_product, c("M_PE-P", "M_Fe3+-cytochrome_b5", "2 M_H2O"))
    
    ## pep_to_lpep
    reaction <- "RHEA:36195"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE-P = M_H+ + M_LPE-P + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36196"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE-P => M_H+ + M_LPE-P + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36197"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE-P <= M_H+ + M_LPE-P + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    reaction <- "RHEA:36198"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PE-P <=> M_H+ + M_LPE-P + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PE-P"))
    expect_equal(template$reaction_product, c("M_H+", "M_LPE-P", "M_FA"))
    
    ## lpep_to_fal
    reaction <- "RHEA:16905"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,"M_H2O + M_1-LPE-P = M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16906"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,"M_H2O + M_1-LPE-P => M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16907"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,"M_H2O + M_1-LPE-P <= M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    reaction <- "RHEA:16908"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula,"M_H2O + M_1-LPE-P <=> M_FAL + M_Glycerophosphoethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_FAL", "M_Glycerophosphoethanolamine"))
    
    ## lpep_to_lpap
    reaction <- "RHEA:36203"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P = M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36204"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P => M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36205"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P <= M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    reaction <- "RHEA:36206"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P <=> M_LPA-P + M_H+ + M_Ethanolamine")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_LPA-P", "M_H+", "M_Ethanolamine"))
    
    ## lpep_to_mgp
    reaction <- "RHEA:36199"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P = M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36200"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P => M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36201"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P <= M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    reaction <- "RHEA:36202"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_1-LPE-P <=> M_Phosphoethanolamine + M_H+ + M_1-MG-P")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_1-LPE-P"))
    expect_equal(template$reaction_product, c("M_Phosphoethanolamine", "M_H+", "M_1-MG-P"))
    
    ## pep_to_napep_sn1
    reaction <- "pep_to_napep_sn1"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-P + M_PC <=> M_NAPEP + M_2-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPEP", "M_2-LPC"))
    
    ## pep_to_napep_sn2
    reaction <- "pep_to_napep_sn2"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_PE-P + M_PC <=> M_NAPEP + M_1-LPC")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_PE-P", "M_PC"))
    expect_equal(template$reaction_product, c("M_NAPEP", "M_1-LPC"))
    
    ## pg_to_cl
    reaction <- "RHEA:32931"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-DG + M_PG = M_H+ + M_CMP + M_CL")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32932"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-DG + M_PG => M_H+ + M_CMP + M_CL")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32933"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-DG + M_PG <= M_H+ + M_CMP + M_CL")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    reaction <- "RHEA:32934"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_CDP-DG + M_PG <=> M_H+ + M_CMP + M_CL")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_CDP-DG", "M_PG"))
    expect_equal(template$reaction_product, c("M_H+", "M_CMP", "M_CL"))
    
    ## pgp_to_pg
    reaction <- "RHEA:33751"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PGP = M_Pi + M_PG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(template$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33752"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PGP => M_Pi + M_PG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(template$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33753"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PGP <= M_Pi + M_PG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(template$reaction_product, c("M_Pi", "M_PG"))
    
    reaction <- "RHEA:33754"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PGP <=> M_Pi + M_PG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PGP"))
    expect_equal(template$reaction_product, c("M_Pi", "M_PG"))
    
    ## pi_to_dg
    reaction <- "RHEA:43484"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI = M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43485"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI => M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43486"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI <= M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    reaction <- "RHEA:43487"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI <=> M_myo-Inositol-1-P + M_H+ + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_myo-Inositol-1-P", "M_H+", "M_1,2-DG"))
    
    ## pi_to_sn1lpi
    reaction <- "RHEA:18001"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI = M_1-LPI + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18002"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI => M_1-LPI + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18003"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI <= M_1-LPI + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    reaction <- "RHEA:18004"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_PI <=> M_1-LPI + M_H+ + M_FA")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_PI"))
    expect_equal(template$reaction_product, c("M_1-LPI", "M_H+", "M_FA"))
    
    ## ps_to_pe
    reaction <- "RHEA:20828"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_PS = M_CO2 + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(template$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20829"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_PS => M_CO2 + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(template$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20830"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_PS <= M_CO2 + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(template$reaction_product, c("M_CO2", "M_PE"))
    
    reaction <- "RHEA:20831"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H+ + M_PS <=> M_CO2 + M_PE")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H+", "M_PS"))
    expect_equal(template$reaction_product, c("M_CO2", "M_PE"))
    
    ## sm_to_cer
    reaction <- "sm_to_cer"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_SM <=> M_Phosphocholine + M_H+ + M_Cer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_SM"))
    expect_equal(template$reaction_product, c("M_Phosphocholine", "M_H+", "M_Cer"))
    
    ## sphinga_to_dhcer
    reaction <- "RHEA:53424"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_Sphinganine = M_CoA + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(template$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53425"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_Sphinganine => M_CoA + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(template$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53426"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_Sphinganine <= M_CoA + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(template$reaction_product, c("M_CoA", "M_DhCer"))
    
    reaction <- "RHEA:53427"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_AcylCoA + M_Sphinganine <=> M_CoA + M_DhCer")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_AcylCoA", "M_Sphinganine"))
    expect_equal(template$reaction_product, c("M_CoA", "M_DhCer"))
    
    ## tg_to_dg
    reaction <- "RHEA:33271"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33272"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33273"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:33274"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG <=> M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44864"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG = M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44865"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG => M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44866"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG <= M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
    
    reaction <- "RHEA:44867"
    template <- .create_template(template = NA, reaction = reaction)
    expect_equal(template$reaction_name, "")
    expect_equal(template$reaction_formula, "M_H2O + M_TG <=> M_H+ + M_FA + M_1,2-DG")
    expect_equal(template$reaction_RHEA, reaction)
    expect_equal(template$reaction_isReversible, "")
    expect_equal(template$reaction_geneAssociation, "")
    expect_equal(template$reaction_pathway, "")
    expect_equal(template$reaction_substrate, c("M_H2O", "M_TG"))
    expect_equal(template$reaction_product, c("M_H+", "M_FA", "M_1,2-DG"))
})

