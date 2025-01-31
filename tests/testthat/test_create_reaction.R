## function .create_reaction
test_that(".create_reaction works", {
    
    ## dg_to_pa
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:10272")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP = PA(18:0/16:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/16:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 = CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:10273")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP => PA(18:0/16:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/16:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 => CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:10274")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP <= PA(18:0/16:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/16:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <= CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:10275")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + ATP <=> PA(18:0/16:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/16:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:30616 <=> CHEBI:58608 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:456216 + CHEBI:15378"))
    
    ## pc_to_dg
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:10604")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(20:0/18:0) + H2O = DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(20:0/18:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(20:0/18:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:10605")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(20:0/18:0) + H2O => DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(20:0/18:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(20:0/18:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:10606")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(20:0/18:0) + H2O <= DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(20:0/18:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(20:0/18:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:10607")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(20:0/18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(20:0/18:0) + H2O <=> DG(20:0/18:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(20:0/18:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(20:0/18:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:15378 + CHEBI:295975"))
    
    ## dg_to_tg
    dg <- "DG(18:0/16:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:10868")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) = TG(18:0/16:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("TG(18:0/16:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 = CHEBI:64615 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64615 + CHEBI:57287"))
    
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) = TG(18:0/16:0/14:0) + CoA")
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:10869")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) => TG(18:0/16:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("TG(18:0/16:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 => CHEBI:64615 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64615 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:10870")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) <= TG(18:0/16:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("TG(18:0/16:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <= CHEBI:64615 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64615 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:10871")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CoA(14:0) <=> TG(18:0/16:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("TG(18:0/16:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58342 <=> CHEBI:64615 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64615 + CHEBI:57287"))
    
    
    ## cdpdg_to_pi
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:11580")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol = PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + myo-Inositol"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 = CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:17268"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:11581")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol => PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + myo-Inositol"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 => CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:17268"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:11582")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol <= PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + myo-Inositol"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <= CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:17268"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:11583")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PI, "PI(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + myo-Inositol <=> PI(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + myo-Inositol"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:17268 <=> CHEBI:57880 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:17268"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:60377 + CHEBI:15378"))
    
    ## cer_to_glccer
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:12088")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$GlcCer, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose = GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose"))
    expect_equal(reaction_l[[2]]$reaction_product, c("GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 = CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:58885"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:22801 + CHEBI:15378 + CHEBI:58223"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:12089")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$GlcCer, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose => GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose"))
    expect_equal(reaction_l[[2]]$reaction_product, c("GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 => CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:58885"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:22801 + CHEBI:15378 + CHEBI:58223"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:12090")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$GlcCer, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose <= GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose"))
    expect_equal(reaction_l[[2]]$reaction_product, c("GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <= CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:58885"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:22801 + CHEBI:15378 + CHEBI:58223"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:12091")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$GlcCer, "GlcCer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose <=> GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + UDP-Glucose"))
    expect_equal(reaction_l[[2]]$reaction_product, c("GlcCer(d16:1(3OH,4OH)(15Me)/18:0) + H+ + UDP"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:58885 <=> CHEBI:22801 + CHEBI:15378 + CHEBI:58223")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:58885"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:22801 + CHEBI:15378 + CHEBI:58223"))
    
    ## cdpdg_to_pgp
    cdpdg <- "CDP-DG(12:0(11Me)/14:0)"
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:12593")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P = PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PGP(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 = CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60110 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:12594")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P => PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PGP(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 => CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60110 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:12595")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P <= PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PGP(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <= CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60110 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CDPDG = cdpdg),
        template = list(), reaction =  "RHEA:12596")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[1]]$PGP, "PGP(12:0(11Me)/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P <=> PGP(12:0(11Me)/14:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CDP-DG(12:0(11Me)/14:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PGP(12:0(11Me)/14:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58332 + CHEBI:57597 <=> CHEBI:60110 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58332 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60110 + CHEBI:60377 + CHEBI:15378"))
    
    ## sn1lpc_to_pc
    sn1lpc <- "PC(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:12937")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) = PC(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 = CHEBI:57643 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:12938")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) => PC(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 => CHEBI:57643 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:12939")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <= PC(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <= CHEBI:57643 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:12940")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + CoA(18:0) <=> PC(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:58342 <=> CHEBI:57643 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:57287"))
    
    ## pc_to_pa
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:14445")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O = PA(16:0/14:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(16:0/14:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:14446")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O => PA(16:0/14:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(16:0/14:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:14447")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <= PA(16:0/14:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(16:0/14:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:14448")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <=> PA(16:0/14:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(16:0/14:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:15354 + CHEBI:15378"))
    
    ## sn1lpc_to_fa
    sn1lpc <- "PC(14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = list(), reaction =  "RHEA:15177")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + H2O = FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = list(), reaction =  "RHEA:15178")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + H2O => FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = list(), reaction =  "RHEA:15179")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + H2O <= FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPC = sn1lpc), 
        template = list(), reaction =  "RHEA:15180")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(14:0/0:0) + H2O <=> FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58168 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58168 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    ## coa_to_lpa
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:15325")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P = PA(18:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 = CHEBI:57970 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:15326")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P => PA(18:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 => CHEBI:57970 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:15327")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P <= PA(18:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <= CHEBI:57970 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:15328")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Glycerol-3-P <=> PA(18:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Glycerol-3-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57597 <=> CHEBI:57970 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57597"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:57287"))
    
    ## fa_to_coa
    fa <- "FA(18:0)"
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:15421")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA = CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 = CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57560 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83139 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:15422")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA => CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 => CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57560 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83139 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:15423")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA <= CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <= CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57560 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83139 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:15424")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA <=> CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57560 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:83139 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57560 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83139 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa),
        template = list(), reaction =  "RHEA:38883")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA = CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 = CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:28868 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77636 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:38884")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA => CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 => CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:28868 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77636 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:38885")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA <= CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <= CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:28868 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77636 + CHEBI:456215 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(FA = fa), 
        template = list(), reaction =  "RHEA:38886")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "FA(18:0) + ATP + CoA <=> CoA(18:0) + AMP + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("FA(18:0) + ATP + CoA"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA(18:0) + AMP + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:28868 + CHEBI:30616 + CHEBI:57287 <=> CHEBI:77636 + CHEBI:456215 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:28868 + CHEBI:30616 + CHEBI:57287"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77636 + CHEBI:456215 + CHEBI:33019"))
    
    ## pc_to_sn1lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:15801")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O = PC(16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:15802")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O => PC(16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:15803")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <= PC(16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:15804")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <=> PC(16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:58168 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:28868 + CHEBI:15378"))
    
    ## pa_to_cdpdg
    pa <- "PA(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:16229")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ = CDP-DG(14:0/16:0) + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + CTP + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CDP-DG(14:0/16:0) + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 = CHEBI:58332 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:37563 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58332 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:16230")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ => CDP-DG(14:0/16:0) + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + CTP + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CDP-DG(14:0/16:0) + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 => CHEBI:58332 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:37563 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58332 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:16231")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ <= CDP-DG(14:0/16:0) + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + CTP + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CDP-DG(14:0/16:0) + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <= CHEBI:58332 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:37563 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58332 + CHEBI:33019"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:16232")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$CDPDG, "CDP-DG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + CTP + H+ <=> CDP-DG(14:0/16:0) + PPi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + CTP + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CDP-DG(14:0/16:0) + PPi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:37563 + CHEBI:15378 <=> CHEBI:58332 + CHEBI:33019")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:37563 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58332 + CHEBI:33019"))
    
    ## lpep_to_pep
    lpep <- "PE(P-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = lpep, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16245")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) = PE(P-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 = CHEBI:77290 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPEP = lpep, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16246")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) => PE(P-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 => CHEBI:77290 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPEP = lpep, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16247")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) <= PE(P-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <= CHEBI:77290 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPEP = lpep, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16248")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + CoA(18:0) <=> PE(P-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:58342 <=> CHEBI:77290 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:57287"))
    
    ## sn2mg_to_dg
    sn2mg <- "MG(0:0/14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16741")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) = DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16742")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) => DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16743")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <= DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:16744")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <=> DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32947")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) = DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32948")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) => DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32949")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <= DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32950")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) + CoA(16:0) <=> DG(16:0/14:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(16:0/14:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17389 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17389 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    ## lpep_to_fal
    sn1lpep <- "PE(P-16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:16905")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FAL, "FAL(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O = FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAL(16:0) + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:73359 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73359 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:16906")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FAL, "FAL(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O => FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAL(16:0) + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:73359 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73359 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:16907")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FAL, "FAL(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <= FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAL(16:0) + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:73359 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73359 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:16908")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FAL, "FAL(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <=> FAL(16:0) + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAL(16:0) + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:73359 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73359 + CHEBI:143890"))
    
    ## nae_to_fa
    nae <- "NAE(18:0)"
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:17505")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) = FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 = CHEBI:57560 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:15897"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57560 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:17506")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) => FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 => CHEBI:57560 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:15897"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57560 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:17507")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) <= FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <= CHEBI:57560 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:15897"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57560 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:17508")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) <=> FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:15897 <=> CHEBI:57560 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:15897"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57560 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:39995")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) = FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 = CHEBI:28868 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:52640"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:39996")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) => FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 => CHEBI:28868 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:52640"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:39997")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) <= FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <= CHEBI:28868 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:52640"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(NAE = nae), 
        template = list(), reaction =  "RHEA:39998")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAE(18:0) <=> FA(18:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAE(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(18:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:52640 <=> CHEBI:28868 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:52640"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:57603"))
    
    ## coa_to_acyldhap
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17657")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P = DHAP(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Dihydroxyacetone-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 = CHEBI:57534 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57642"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57534 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17658")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P => DHAP(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Dihydroxyacetone-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 => CHEBI:57534 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57642"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57534 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17659")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P <= DHAP(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Dihydroxyacetone-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <= CHEBI:57534 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57642"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57534 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17660")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + Dihydroxyacetone-P <=> DHAP(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + Dihydroxyacetone-P"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 + CHEBI:57642 <=> CHEBI:57534 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + CHEBI:57642"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57534 + CHEBI:57287"))
    
    ## coa_to_ce
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17729")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$CE, "CE(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cholesterol + CoA(18:0) = CE(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c( "Cholesterol + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CE(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 = CHEBI:17002 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:16113 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17002 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17730")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$CE, "CE(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cholesterol + CoA(18:0) => CE(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c( "Cholesterol + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CE(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 => CHEBI:17002 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:16113 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17002 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17731")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$CE, "CE(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cholesterol + CoA(18:0) <= CE(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c( "Cholesterol + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CE(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 <= CHEBI:17002 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:16113 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17002 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:17732")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$CE, "CE(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cholesterol + CoA(18:0) <=> CE(18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c( "Cholesterol + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CE(18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:16113 + CHEBI:58342 <=> CHEBI:17002 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:16113 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17002 + CHEBI:57287"))
    
    ## cer_to_cerp
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:17929")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP = ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 = CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:456216 + CHEBI:57674 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:17930")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP => ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 => CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:456216 + CHEBI:57674 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:17931")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP <= ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <= CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:456216 + CHEBI:57674 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:17932")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP <=> ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:1(3OH,4OH)(15Me)/18:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("ADP + CerP(d16:1(3OH,4OH)(15Me)/18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52639 + CHEBI:30616 <=> CHEBI:456216 + CHEBI:57674 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52639 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:456216 + CHEBI:57674 + CHEBI:15378"))
    
    ## pi_to_sn1lpi
    pi <- "PI(16:0/18:1(9Z))"
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:18001")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O = PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/0:0) + FA(18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64771 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:18002")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O => PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/0:0) + FA(18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64771 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:18003")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O <= PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/0:0) + FA(18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64771 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:18004")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O <=> PI(16:0/0:0) + FA(18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/0:0) + FA(18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:64771 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64771 + CHEBI:28868 + CHEBI:15378"))
    
    ## pc_to_sn2lpc
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:18689")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O = PC(0:0/14:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/14:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 = CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:18690")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O => PC(0:0/14:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/14:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 => CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:18691")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <= PC(0:0/14:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/14:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <= CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:18692")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + H2O <=> PC(0:0/14:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/14:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:15377 <=> CHEBI:57875 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:28868 + CHEBI:15378"))
    
    ## cer_to_sm
    cer <- "Cer(d16:1(3OH,4OH)(15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:18765")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) = DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 = CHEBI:17815 + CHEBI:17636")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:52639"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:17636"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:18766")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) => DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 => CHEBI:17815 + CHEBI:17636")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:52639"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:17636"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:18767")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) <= DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <= CHEBI:17815 + CHEBI:17636")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:52639"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:17636"))
    
    reaction_l <- .create_reaction(substrates = list(Cer = cer),
        template = list(), reaction =  "RHEA:18768")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "PC + Cer(d16:1(3OH,4OH)(15Me)/18:0) <=> DG + SM(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:52639 <=> CHEBI:17815 + CHEBI:17636")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:52639"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:17636"))
    
    ## lpa_to_pa
    lpa <- "PA(18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(
        substrates = list(sn1LPA = lpa, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:19709")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) = PA(18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 = CHEBI:58608 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57970 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPA = lpa, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:19710")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) => PA(18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 => CHEBI:58608 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57970 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPA = lpa, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:19711")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <= PA(18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <= CHEBI:58608 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57970 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPA = lpa, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:19712")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(18:0/0:0) + CoA(14:0) <=> PA(18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57970 + CHEBI:58342 <=> CHEBI:58608 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57970 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:57287"))
    
    ## ps_to_pe
    ps <- "PS(14:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:20828")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PS(14:0/14:0) + H+ = PE(14:0/14:0) + CO2")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/14:0) + CO2"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 = CHEBI:64612 + CHEBI:16526")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:16526"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:20829")
    ps <- "PS(14:0/14:0)"
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PS(14:0/14:0) + H+ => PE(14:0/14:0) + CO2")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/14:0) + CO2"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 => CHEBI:64612 + CHEBI:16526")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:16526"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:20830")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PS(14:0/14:0) + H+ <= PE(14:0/14:0) + CO2")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/14:0) + CO2"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <= CHEBI:64612 + CHEBI:16526")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:16526"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:20831")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/14:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PS(14:0/14:0) + H+ <=> PE(14:0/14:0) + CO2")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/14:0) + CO2"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15378 <=> CHEBI:64612 + CHEBI:16526")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:16526"))
    
    ## peo_to_pep
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = list(), reaction =  "RHEA:22956")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 = PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377"))
    
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = list(), reaction =  "RHEA:22957")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 => PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377"))
    
    
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = list(), reaction =  "RHEA:22958")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <= PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377"))
    
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = list(), reaction =  "RHEA:22959")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2 <=> PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/14:0) + Fe3+-cytochrome b5 + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77290 + CHEBI:29034 + 2 CHEBI:15377"))
    
    ## lpco_to_pco
    sn1lpco <- "PC(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(
        substrates = list(sn1LPCO = sn1lpco, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:23992")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) = PC(O-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 = CHEBI:36702 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPCO = sn1lpco, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:23993")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) => PC(O-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 => CHEBI:36702 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPCO = sn1lpco, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:23994")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <= PC(O-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <= CHEBI:36702 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPCO = sn1lpco, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:23995")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + CoA(18:0) <=> PC(O-16:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:58342 <=> CHEBI:36702 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:57287"))
    
    ## pa_to_dg
    pa <- "PA(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:27429")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + H2O = DG(14:0/16:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 = CHEBI:17815 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:27430")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + H2O => DG(14:0/16:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 => CHEBI:17815 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:27431")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + H2O <= DG(14:0/16:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <= CHEBI:17815 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PA = pa), 
        template = list(), reaction = "RHEA:27432")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(14:0/16:0) + H2O <=> DG(14:0/16:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58608 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58608 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:43474"))
    
    ## pe_to_ps
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:27606")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + L-Serine = PS(14:0/16:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(14:0/16:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 = CHEBI:57262 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:27607")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + L-Serine => PS(14:0/16:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(14:0/16:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 => CHEBI:57262 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:27608")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + L-Serine <= PS(14:0/16:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(14:0/16:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <= CHEBI:57262 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:57603"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:27609")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + L-Serine <=> PS(14:0/16:0) + Ethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(14:0/16:0) + Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:57603")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:57603"))

    ## sn2mg_to_fa
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = list(), reaction =  "RHEA:32871")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) = FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + MG(0:0/14:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 = CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:17389"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:29067 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = list(), reaction =  "RHEA:32872")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) => FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + MG(0:0/14:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 => CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:17389"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:29067 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = list(), reaction =  "RHEA:32873")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <= FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + MG(0:0/14:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <= CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:17389"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:29067 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = list(), reaction =  "RHEA:32874")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + MG(0:0/14:0/0:0) <=> FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + MG(0:0/14:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:17389 <=> CHEBI:29067 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:17389"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:29067 + CHEBI:17754 + CHEBI:15378"))
    
    ## pg_to_cl
    pg <- "PG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    cdpdg <- "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = list(), reaction = "RHEA:32931")
    expect_equal(reaction_l[[1]]$PG, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CDPDG, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 = CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:58332"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = list(), reaction = "RHEA:32932")
    expect_equal(reaction_l[[1]]$PG, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CDPDG, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 => CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:58332"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = list(), reaction = "RHEA:32933")
    expect_equal(reaction_l[[1]]$PG, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CDPDG, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <= CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:58332"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg, CDPDG = cdpdg), 
        template = list(), reaction = "RHEA:32934")
    expect_equal(reaction_l[[1]]$PG, 
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CDPDG, 
        "CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)")
    expect_equal(reaction_l[[1]]$CL, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0])")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(18:4(6Z,9Z,12Z,15Z)/14:0) + CDP-DG(18:4(6Z,9Z,12Z,15Z)/14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/14:0],3'-[18:4(6Z,9Z,12Z,15Z)/14:0]) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:58332 <=> CHEBI:62237 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:58332"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:60377 + CHEBI:15378"))
    
    ## cl_to_lcl
    cl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])"
    reaction_l <- .create_reaction(substrates = list(CL = cl),
        template = list(), reaction = "RHEA:32935")
    expect_equal(reaction_l[[1]]$CL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$LCL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 = CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:62237 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64743 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CL = cl),
        template = list(), reaction = "RHEA:32936")
    expect_equal(reaction_l[[1]]$CL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$LCL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 => CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:62237 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64743 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CL = cl),
        template = list(), reaction = "RHEA:32937")
    expect_equal(reaction_l[[1]]$CL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$LCL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <= CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:62237 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64743 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(CL = cl),
        template = list(), reaction = "RHEA:32938")
    expect_equal(reaction_l[[1]]$CL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$LCL,
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$FA, "FA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + FA(18:4(6Z,9Z,12Z,15Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:62237 + CHEBI:15377 <=> CHEBI:64743 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:62237 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64743 + CHEBI:28868 + CHEBI:15378"))
    
    ## dg_to_pc
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32939")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline = PC(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 = CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32940")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline => PC(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 => CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32941")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline <= PC(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <= CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32942")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Choline <=> PC(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:58779 <=> CHEBI:57643 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57643 + CHEBI:60377 + CHEBI:15378"))
    
    ## dg_to_pe
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32943")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine = PE(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 = CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32944")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine => PE(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 => CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:32945")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine <= PE(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <= CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg), 
        template = list(), reaction =  "RHEA:32946")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + CDP-Ethanolamine <=> PE(18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:57876 <=> CHEBI:64612 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:60377 + CHEBI:15378"))
    
    ## sn1lpe_to_fa
    pe <- "PE(14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = list(), reaction =  "RHEA:32967")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + H2O = FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = list(), reaction =  "RHEA:32968")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + H2O => FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe), 
        template = list(), reaction =  "RHEA:32969")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + H2O <= FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:143890"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = pe),
        template = list(), reaction =  "RHEA:32970")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + H2O <=> FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:143890"))
    
    ## sn1lpe_to_pe
    sn1lpe <- "PE(14:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32995")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) = PE(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 = CHEBI:64612 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32996")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) => PE(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 => CHEBI:64612 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32997")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) <= PE(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <= CHEBI:64612 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPE = sn1lpe, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:32998")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/0:0) + CoA(18:0) <=> PE(14:0/18:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/0:0) + CoA(18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/18:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64381 + CHEBI:58342 <=> CHEBI:64612 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64381 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64612 + CHEBI:57287"))
    
    ## nape_to_nae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:33159")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) = PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/16:0) + NAE(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:52640 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:33160")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) => PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/16:0) + NAE(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:52640 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:33161")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <= PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/16:0) + NAE(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:52640 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:33162")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$NAE, "NAE(18:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <=> PA(14:0/16:0) + NAE(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/16:0) + NAE(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:58608 + CHEBI:52640 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:52640 + CHEBI:15378"))
    
    ## sn1lpi_to_pi
    sn1lpi <- "PI(16:0/0:0)"
    acylcoa <- "CoA(18:1(9Z))"
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:33195")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) = PI(16:0/18:1(9Z)) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/0:0) + CoA(18:1(9Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/18:1(9Z)) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 = CHEBI:57880 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64771 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:33196")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) => PI(16:0/18:1(9Z)) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/0:0) + CoA(18:1(9Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/18:1(9Z)) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 => CHEBI:57880 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64771 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:33197")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <= PI(16:0/18:1(9Z)) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/0:0) + CoA(18:1(9Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/18:1(9Z)) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <= CHEBI:57880 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64771 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(sn1LPI = sn1lpi, AcylCoA = acylcoa), 
        template = list(), reaction = "RHEA:33198")
    expect_equal(reaction_l[[1]]$sn1LPI, "PI(16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:1(9Z))")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PI(16:0/0:0) + CoA(18:1(9Z)) <=> PI(16:0/18:1(9Z)) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/0:0) + CoA(18:1(9Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PI(16:0/18:1(9Z)) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64771 + CHEBI:58342 <=> CHEBI:57880 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64771 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57880 + CHEBI:57287"))
    
    ## sn1lpg_to_pg
    sn1lpg <- "PG(14:0/0:0)"
    acylcoa <- "CoA(16:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPG = sn1lpg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:33203")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/0:0) + CoA(16:0) = PG(14:0/16:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/16:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 = CHEBI:64716 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64840 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPG = sn1lpg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:33204")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/0:0) + CoA(16:0) => PG(14:0/16:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/16:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 => CHEBI:64716 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64840 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPG = sn1lpg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:33205")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/0:0) + CoA(16:0) <= PG(14:0/16:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/16:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 <= CHEBI:64716 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64840 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPG = sn1lpg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:33206")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/0:0) + CoA(16:0) <=> PG(14:0/16:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/16:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64840 + CHEBI:58342 <=> CHEBI:64716 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64840 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:57287"))
    
    ## tg_to_dg
    tg <- "TG(18:0/16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:33271")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O = DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O = DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 = CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 = CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:33272")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O => DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O => DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 => CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 => CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:33273")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O <= DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O <= DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 <= CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 <= CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:33274")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O <=> DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O <=> DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 <=> CHEBI:17815 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:17815 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:44864")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O = DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O = DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 = CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 = CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:44865")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O => DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O => DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 => CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 => CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:44866")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O <= DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O <= DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 <= CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 <= CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(TG = tg), 
        template = list(), reaction = "RHEA:44867")
    expect_equal(reaction_l[[1]]$TG, "TG(18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1Loss_FA, "FA(18:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn3Loss_FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula[1],
        "TG(18:0/16:0/14:0) + H2O <=> DG(14:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_formula[2],
        "TG(18:0/16:0/14:0) + H2O <=> DG(18:0/16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible[1], "")
    expect_equal(reaction_l[[2]]$reaction_isReversible[2], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[1], "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation[2], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[1], "")
    expect_equal(reaction_l[[2]]$reaction_pathway[2], "")
    expect_equal(reaction_l[[2]]$reaction_substrate[1], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_substrate[2], c("TG(18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product[1], c("DG(14:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_product[2], c("DG(18:0/16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[1], "CHEBI:64615 + CHEBI:15377 <=> CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_formula_chebi[2], "CHEBI:64615 + CHEBI:15377 <=> CHEBI:49172 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[1], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi[2], c("CHEBI:64615 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[1], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi[2], c("CHEBI:49172 + CHEBI:28868 + CHEBI:15378"))
    
    ## dg_to_sn2mg
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:33275")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(0:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 = CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17389 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:33276")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(0:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 => CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17389 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:33277")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(0:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <= CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17389 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:33278")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(0:0/16:0/0:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(0:0/16:0/0:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <=> CHEBI:17389 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17389 + CHEBI:28868 + CHEBI:15378"))
    
    ## cerp_to_cer
    cerp <- "CerP(d16:1(3OH,4OH)(15Me)/18:0)"
    reaction_l <- .create_reaction(substrates = list(CerP = cerp),
        template = list(), reaction =  "RHEA:33743")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) = Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 = CHEBI:43474 + CHEBI:52639")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:57674"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:43474 + CHEBI:52639"))
    
    reaction_l <- .create_reaction(substrates = list(CerP = cerp),
        template = list(), reaction =  "RHEA:33744")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) => Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 => CHEBI:43474 + CHEBI:52639")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:57674"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:43474 + CHEBI:52639"))
    
    reaction_l <- .create_reaction(substrates = list(CerP = cerp),
        template = list(), reaction =  "RHEA:33745")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) <= Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <= CHEBI:43474 + CHEBI:52639")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:57674"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:43474 + CHEBI:52639"))
    
    reaction_l <- .create_reaction(substrates = list(CerP = cerp),
        template = list(), reaction =  "RHEA:33746")
    expect_equal(reaction_l[[1]]$CerP, "CerP(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0) <=> Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + CerP(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Pi + Cer(d16:1(3OH,4OH)(15Me)/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:57674 <=> CHEBI:43474 + CHEBI:52639")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:57674"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:43474 + CHEBI:52639"))
    
    ## sn1mg_to_lpa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:33747")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP = PA(14:0/0:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/0:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 = CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:33748")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP => PA(14:0/0:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/0:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 => CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:33749")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP <= PA(14:0/0:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/0:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <= CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:456216 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:33750")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPA, "PA(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + ATP <=> PA(14:0/0:0) + ADP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + ATP"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/0:0) + ADP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:30616 <=> CHEBI:57970 + CHEBI:456216 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:30616"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57970 + CHEBI:456216 + CHEBI:15378"))
    
    ## pgp_to_pg
    pgp <- "PGP(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = list(), reaction = "RHEA:33751")
    expect_equal(reaction_l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PGP(16:0/14:0) + H2O = PG(16:0/14:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PGP(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(16:0/14:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 = CHEBI:64716 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:60110 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = list(), reaction = "RHEA:33752")
    expect_equal(reaction_l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PGP(16:0/14:0) + H2O => PG(16:0/14:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PGP(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(16:0/14:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 => CHEBI:64716 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:60110 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = list(), reaction = "RHEA:33753")
    expect_equal(reaction_l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PGP(16:0/14:0) + H2O <= PG(16:0/14:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PGP(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(16:0/14:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <= CHEBI:64716 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:60110 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PGP = pgp), 
        template = list(), reaction = "RHEA:33754")
    expect_equal(reaction_l[[1]]$PGP, "PGP(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PG, "PG(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PGP(16:0/14:0) + H2O <=> PG(16:0/14:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PGP(16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(16:0/14:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:60110 + CHEBI:15377 <=> CHEBI:64716 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:60110 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64716 + CHEBI:43474"))
    
    ## sn1mg_to_fa
    sn1mg <- "MG(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:34019")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O = FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 = CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:34020")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O => FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 => CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:34021")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O <= FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <= CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:17754 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg), 
        template = list(), reaction =  "RHEA:34022")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + H2O <=> FA(14:0) + Glycerol + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + Glycerol + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:17754 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:17754 + CHEBI:15378"))
    
    ## dg_to_sn1mg
    dg <- "DG(18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:35663")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 = CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64683 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:35664")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 => CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64683 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:35665")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <= CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64683 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:35666")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:17815 + CHEBI:15377 <=> CHEBI:64683 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:17815 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64683 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:44712")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O = MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 = CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:49172 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:35759 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:44713")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O => MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 => CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:49172 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:35759 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:44714")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <= MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <= CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:49172 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:35759 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DG = dg),
        template = list(), reaction =  "RHEA:44715")
    expect_equal(reaction_l[[1]]$DG, "DG(18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(18:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(18:0/16:0/0:0) + H2O <=> MG(18:0/0:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(18:0/16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(18:0/0:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:49172 + CHEBI:15377 <=> CHEBI:35759 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:49172 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:35759 + CHEBI:28868 + CHEBI:15378"))

    ## lcl_to_cl
    lcl <- "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])"
    acylcoa <- "CoA(18:4(6Z,9Z,12Z,15Z))"
    reaction_l <- .create_reaction(
        substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:35839")
    expect_equal(reaction_l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) = CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 = CHEBI:62237 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64743 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:35840")
    expect_equal(reaction_l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) => CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 => CHEBI:62237 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64743 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:35841")
    expect_equal(reaction_l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <= CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <= CHEBI:62237 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64743 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(
        substrates = list(LCL = lcl, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:35842")
    expect_equal(reaction_l[[1]]$LCL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:4(6Z,9Z,12Z,15Z))")
    expect_equal(reaction_l[[1]]$CL, "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)])")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z)) <=> CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[0:0/18:4(6Z,9Z,12Z,15Z)]) + CoA(18:4(6Z,9Z,12Z,15Z))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CL(1'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)],3'-[18:4(6Z,9Z,12Z,15Z)/18:4(6Z,9Z,12Z,15Z)]) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64743 + CHEBI:58342 <=> CHEBI:62237 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64743 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:62237 + CHEBI:57287"))
    
    ## lpco_to_mgo
    sn1lpco <- "PC(O-16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:36083")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O = MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(O-16:0/0:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15850 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:36084")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O => MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(O-16:0/0:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15850 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:36085")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <= MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(O-16:0/0:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15850 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:36086")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGO, "MG(O-16:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <=> MG(O-16:0/0:0/0:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(O-16:0/0:0/0:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:15850 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15850 + CHEBI:15378 + CHEBI:295975"))
    
    ## acyldhap_to_alkyldhap
    acyldhap <- "DHAP(18:0)"
    fao <- "FAO(16:0)"
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FAO = fao),
        template = list(), reaction = "RHEA:36171")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(16:0)")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(18:0) + FAO(16:0) = DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(18:0) + FAO(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(O-16:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 = CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57534 + CHEBI:17135"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73315 + CHEBI:57560 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FAO = fao),
        template = list(), reaction = "RHEA:36172")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(16:0)")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(18:0) + FAO(16:0) => DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(18:0) + FAO(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(O-16:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 => CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57534 + CHEBI:17135"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73315 + CHEBI:57560 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FAO = fao),
        template = list(), reaction = "RHEA:36173")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(16:0)")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(18:0) + FAO(16:0) <= DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(18:0) + FAO(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(O-16:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <= CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57534 + CHEBI:17135"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73315 + CHEBI:57560 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylDHAP = acyldhap, FAO = fao),
        template = list(), reaction = "RHEA:36174")
    expect_equal(reaction_l[[1]]$AcylDHAP, "DHAP(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(16:0)")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(18:0) + FAO(16:0) <=> DHAP(O-16:0) + FA(18:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(18:0) + FAO(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DHAP(O-16:0) + FA(18:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57534 + CHEBI:17135 <=> CHEBI:73315 + CHEBI:57560 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57534 + CHEBI:17135"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73315 + CHEBI:57560 + CHEBI:15378"))
    
    ## alkyldhap_to_lpao
    alkyldhap <- "DHAP(O-18:0)"
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap),
        template = list(), reaction =  "RHEA:36175")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(O-18:0) + H+ + NADPH = PA(O-18:0/0:0) + NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(O-18:0) + H+ + NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/0:0) + NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 = CHEBI:58014 + CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73315 + CHEBI:15378 + CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap),
        template = list(), reaction =  "RHEA:36176")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(O-18:0) + H+ + NADPH => PA(O-18:0/0:0) + NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(O-18:0) + H+ + NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/0:0) + NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 => CHEBI:58014 + CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73315 + CHEBI:15378 + CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap),
        template = list(), reaction =  "RHEA:36177")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(O-18:0) + H+ + NADPH <= PA(O-18:0/0:0) + NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(O-18:0) + H+ + NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/0:0) + NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <= CHEBI:58014 + CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73315 + CHEBI:15378 + CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AlkylDHAP = alkyldhap),
        template = list(), reaction =  "RHEA:36178")
    expect_equal(reaction_l[[1]]$AlkylDHAP, "DHAP(O-18:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, "DHAP(O-18:0) + H+ + NADPH <=> PA(O-18:0/0:0) + NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DHAP(O-18:0) + H+ + NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/0:0) + NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73315 + CHEBI:15378 + CHEBI:57783 <=> CHEBI:58014 + CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73315 + CHEBI:15378 + CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:58349"))
    
    ## dgo_to_pco
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36179")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "DG(O-18:0/16:0/0:0) + CDP-Choline = PC(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 = CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36180")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "DG(O-18:0/16:0/0:0) + CDP-Choline => PC(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 => CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36181")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "DG(O-18:0/16:0/0:0) + CDP-Choline <= PC(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <= CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36182")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "DG(O-18:0/16:0/0:0) + CDP-Choline <=> PC(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Choline"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:58779 <=> CHEBI:36702 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:58779"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:36702 + CHEBI:60377 + CHEBI:15378"))
    
    ## dgo_to_peo
    dgo <- "DG(O-18:0/16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36187")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine = PE(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 = CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60520 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36188")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine => PE(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 => CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60520 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36189")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine <= PE(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <= CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60520 + CHEBI:60377 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(DGO = dgo), 
        template = list(), reaction =  "RHEA:36190")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-18:0/16:0/0:0)")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "DG(O-18:0/16:0/0:0) + CDP-Ethanolamine <=> PE(O-18:0/16:0) + CMP + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("DG(O-18:0/16:0/0:0) + CDP-Ethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(O-18:0/16:0) + CMP + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:52595 + CHEBI:57876 <=> CHEBI:60520 + CHEBI:60377 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:52595 + CHEBI:57876"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:60520 + CHEBI:60377 + CHEBI:15378"))
    
    ## pep_to_lpep
    pep <- "PE(P-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = list(), reaction =  "RHEA:36195")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + H2O = PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi,  "CHEBI:77290 + CHEBI:15377 = CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77288 + CHEBI:29067 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = list(), reaction =  "RHEA:36196")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + H2O => PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi,  "CHEBI:77290 + CHEBI:15377 => CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77288 + CHEBI:29067 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = list(), reaction =  "RHEA:36197")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + H2O <= PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi,  "CHEBI:77290 + CHEBI:15377 <= CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77288 + CHEBI:29067 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep), 
        template = list(), reaction =  "RHEA:36198")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/14:0) + H2O <=> PE(P-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(P-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi,  "CHEBI:77290 + CHEBI:15377 <=> CHEBI:77288 + CHEBI:29067 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77288 + CHEBI:29067 + CHEBI:15378"))
    
    ## lpep_to_mgp
    sn1lpep <- "PE(P-18:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction = "RHEA:36199")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-18:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGP, "MG(P-18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-18:0/0:0) + H2O = MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-18:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77297 + CHEBI:15378 + CHEBI:58190"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction = "RHEA:36200")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-18:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGP, "MG(P-18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-18:0/0:0) + H2O => MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-18:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77297 + CHEBI:15378 + CHEBI:58190"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction = "RHEA:36201")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-18:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGP, "MG(P-18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-18:0/0:0) + H2O <= MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-18:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77297 + CHEBI:15378 + CHEBI:58190"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction = "RHEA:36202")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-18:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MGP, "MG(P-18:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-18:0/0:0) + H2O <=> MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-18:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(P-18:0/0:0/0:0) + H+ + Phosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77297 + CHEBI:15378 + CHEBI:58190")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77297 + CHEBI:15378 + CHEBI:58190"))
    
    ## lpep_to_lpap
    sn1lpep <- "PE(P-16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:36203")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O = PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(P-16:0/0:0) + Ethanolamine + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 = CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77283 + CHEBI:57603 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:36204")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O => PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(P-16:0/0:0) + Ethanolamine + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 => CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77283 + CHEBI:57603 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:36205")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <= PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(P-16:0/0:0) + Ethanolamine + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <= CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77283 + CHEBI:57603 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPEP = sn1lpep), 
        template = list(), reaction =  "RHEA:36206")
    expect_equal(reaction_l[[1]]$sn1LPEP, "PE(P-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAP, "PA(P-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(P-16:0/0:0) + H2O <=> PA(P-16:0/0:0) + Ethanolamine + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(P-16:0/0:0) + Ethanolamine + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77288 + CHEBI:15377 <=> CHEBI:77283 + CHEBI:57603 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77288 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77283 + CHEBI:57603 + CHEBI:15378"))
    
    ## pco_to_lpco
    pco <- "PC(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = list(), reaction = "RHEA:36231")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/14:0) + H2O = PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 = CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:36702 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:30909 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = list(), reaction = "RHEA:36232")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/14:0) + H2O => PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 => CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:36702 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:30909 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = list(), reaction = "RHEA:36233")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/14:0) + H2O <= PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <= CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:36702 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:30909 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PCO = pco), 
        template = list(), reaction = "RHEA:36234")
    expect_equal(reaction_l[[1]]$PCO, "PC(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/14:0) + H2O <=> PC(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(O-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:36702 + CHEBI:15377 <=> CHEBI:30909 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:36702 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:30909 + CHEBI:28868 + CHEBI:15378"))
    
    ## lpao_to_pao
    lpao <- "PA(O-18:0/0:0)"
    acylcoa <- "CoA(14:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPAO = lpao, AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:36235")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) = PA(O-18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 = CHEBI:73332 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58014 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73332 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPAO = lpao, AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:36236")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) => PA(O-18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 => CHEBI:73332 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58014 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73332 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPAO = lpao,  AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:36237")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <= PA(O-18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <= CHEBI:73332 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58014 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73332 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPAO = lpao,  AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:36238")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-18:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-18:0/0:0) + CoA(14:0) <=> PA(O-18:0/14:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-18:0/0:0) + CoA(14:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-18:0/14:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58014 + CHEBI:58342 <=> CHEBI:73332 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58014 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:73332 + CHEBI:57287"))
    
    ## pao_to_dgo
    pao <- "PA(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = list(), reaction = "RHEA:36239")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-16:0/14:0) + H2O = DG(O-16:0/14:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(O-16:0/14:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 = CHEBI:52595 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73332 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:52595 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = list(), reaction = "RHEA:36240")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-16:0/14:0) + H2O => DG(O-16:0/14:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(O-16:0/14:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 => CHEBI:52595 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73332 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:52595 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = list(), reaction = "RHEA:36241")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-16:0/14:0) + H2O <= DG(O-16:0/14:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(O-16:0/14:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <= CHEBI:52595 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73332 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:52595 + CHEBI:43474"))
    
    reaction_l <- .create_reaction(substrates = list(PAO = pao), 
        template = list(), reaction = "RHEA:36242")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$DGO, "DG(O-16:0/14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PA(O-16:0/14:0) + H2O <=> DG(O-16:0/14:0/0:0) + Pi")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PA(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(O-16:0/14:0/0:0) + Pi"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:73332 + CHEBI:15377 <=> CHEBI:52595 + CHEBI:43474")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:73332 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:52595 + CHEBI:43474"))
    
    ## sn1mg_to_dg
    sn1mg <- "MG(14:0/0:0/0:0)"
    acylcoa <- "CoA(16:0)"
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:38463")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 = CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:38464")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 => CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:38465")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <= CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:38466")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64683 + CHEBI:58342 <=> CHEBI:17815 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64683 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:39943")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) = DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 = CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:39944")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) => DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 => CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))

    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:39945")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <= DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <= CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    reaction_l <- .create_reaction(substrates = list(sn1MG = sn1mg, AcylCoA = acylcoa), 
        template = list(), reaction =  "RHEA:39946")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(14:0/0:0/0:0) + CoA(16:0) <=> DG(14:0/16:0/0:0) + CoA")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(14:0/0:0/0:0) + CoA(16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG(14:0/16:0/0:0) + CoA"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 + CHEBI:58342 <=> CHEBI:49172 + CHEBI:57287")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759 + CHEBI:58342"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:49172 + CHEBI:57287"))
    
    ## lpco_to_lpao
    sn1lpco <- "PC(O-16:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:39927")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O = PA(O-16:0/0:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-16:0/0:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 = CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:39928")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O => PA(O-16:0/0:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-16:0/0:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 => CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:39929")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <= PA(O-16:0/0:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-16:0/0:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <= CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:15354 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(sn1LPCO = sn1lpco), 
        template = list(), reaction =  "RHEA:39930")
    expect_equal(reaction_l[[1]]$sn1LPCO, "PC(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1LPAO, "PA(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(O-16:0/0:0) + H2O <=> PA(O-16:0/0:0) + Choline + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(O-16:0/0:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(O-16:0/0:0) + Choline + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:30909 + CHEBI:15377 <=> CHEBI:58014 + CHEBI:15354 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:30909 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58014 + CHEBI:15354 + CHEBI:15378"))
    
    ## ps_to_sn2lps
    ps <- "PS(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:42212")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPS, "PS(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PS(14:0/16:0) + H2O = PS(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15377 = CHEBI:65214 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65214 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:42213")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPS, "PS(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PS(14:0/16:0) + H2O => PS(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15377 => CHEBI:65214 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65214 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:42214")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPS, "PS(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PS(14:0/16:0) + H2O <= PS(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15377 <= CHEBI:65214 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65214 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PS = ps), 
        template = list(), reaction = "RHEA:42215")
    expect_equal(reaction_l[[1]]$PS, "PS(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPS, "PS(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PS(14:0/16:0) + H2O <=> PS(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PS(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57262 + CHEBI:15377 <=> CHEBI:65214 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57262 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65214 + CHEBI:28868 + CHEBI:15378"))
    
    ## pi_to_dg
    pi <- "PI(16:0/18:1(9Z))"
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:43484")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O = myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 = CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58433 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:43485")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O => myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 => CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58433 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:43486")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O <= myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <= CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58433 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PI = pi), 
        template = list(), reaction = "RHEA:43487")
    expect_equal(reaction_l[[1]]$PI, "PI(16:0/18:1(9Z))")
    expect_equal(reaction_l[[1]]$DG, "DG(16:0/18:1(9Z))")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PI(16:0/18:1(9Z)) + H2O <=> myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PI(16:0/18:1(9Z)) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("myo-Inositol-1-P + DG(16:0/18:1(9Z)) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57880 + CHEBI:15377 <=> CHEBI:58433 + CHEBI:17815 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57880 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58433 + CHEBI:17815 + CHEBI:15378"))
    
    ## pe_to_sn2lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44408")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O = PE(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65213 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44409")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O => PE(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65213 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44410")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O <= PE(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65213 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44411")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/16:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PE(0:0/16:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(0:0/16:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:65213 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:65213 + CHEBI:28868 + CHEBI:15378"))
    
    ## pg_to_sn1lpe
    pg <- "PG(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PG = pg), 
        template = list(), reaction = "RHEA:44416")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/16:0) + H2O = PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 = CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64840 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg), 
        template = list(), reaction = "RHEA:44417")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/16:0) + H2O => PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 => CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64840 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg), 
        template = list(), reaction = "RHEA:44418")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/16:0) + H2O <= PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 <= CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64840 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PG = pg), 
        template = list(), reaction = "RHEA:44419")
    expect_equal(reaction_l[[1]]$PG, "PG(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPG, "PG(14:0/0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PG(14:0/16:0) + H2O <=> PG(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PG(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PG(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64716 + CHEBI:15377 <=> CHEBI:64840 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64716 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64840 + CHEBI:28868 + CHEBI:15378"))
    
    ## pe_to_sn1lpe
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44604")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O = PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 = CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64381 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44605")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O => PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 => CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64381 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44606")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O <= PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <= CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64381 + CHEBI:28868 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:44607")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$sn1LPE, "PE(14:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PE(14:0/0:0) + FA(16:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(14:0/0:0) + FA(16:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:64381 + CHEBI:28868 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:64381 + CHEBI:28868 + CHEBI:15378"))
    
    ## dhcer_to_dhsm
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:44620")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) = DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 = CHEBI:17815 + CHEBI:67090")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:31488"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:67090"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:44621")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) => DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 => CHEBI:17815 + CHEBI:67090")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:31488"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:67090"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:44622")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) <= DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <= CHEBI:17815 + CHEBI:67090")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:31488"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:67090"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:44623")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC + Cer(d16:0(3OH,4OH)(15Me)/12:0) <=> DG + SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC + Cer(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("DG + SM(d16:0(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:31488 <=> CHEBI:17815 + CHEBI:67090")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:31488"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17815 + CHEBI:67090"))
    
    ## sn2lpc_to_fa
    sn2lpc <- "PC(0:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = list(), reaction =  "RHEA:44696")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(0:0/14:0) + H2O = FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(0:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 = CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57875 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = list(), reaction =  "RHEA:44697")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(0:0/14:0) + H2O => FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(0:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 => CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57875 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = list(), reaction =  "RHEA:44698")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(0:0/14:0) + H2O <= FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(0:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <= CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57875 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    reaction_l <- .create_reaction(substrates = list(sn2LPC = sn2lpc), 
        template = list(), reaction =  "RHEA:44699")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(0:0/14:0) + H2O <=> FA(14:0) + H+ + Glycerophosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(0:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57875 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:16870")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57875 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:16870"))
    
    ## pc_to_ps
    pc <- "PC(16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:45088")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + L-Serine = PS(16:0/14:0) + Choline")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(16:0/14:0) + Choline"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 = CHEBI:57262 + CHEBI:15354")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:15354"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:45089")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + L-Serine => PS(16:0/14:0) + Choline")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(16:0/14:0) + Choline"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 => CHEBI:57262 + CHEBI:15354")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:15354"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:45090")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + L-Serine <= PS(16:0/14:0) + Choline")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(16:0/14:0) + Choline"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <= CHEBI:57262 + CHEBI:15354")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:15354"))
    
    reaction_l <- .create_reaction(substrates = list(PC = pc), 
        template = list(), reaction = "RHEA:45091")
    expect_equal(reaction_l[[1]]$PC, "PC(16:0/14:0)")
    expect_equal(reaction_l[[1]]$PS, "PS(16:0/14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(16:0/14:0) + L-Serine <=> PS(16:0/14:0) + Choline")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(16:0/14:0) + L-Serine"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PS(16:0/14:0) + Choline"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:33384 <=> CHEBI:57262 + CHEBI:15354")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:33384"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57262 + CHEBI:15354"))
    
    ## pe_to_nape_sn1
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45188")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) = PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45189")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) => PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45190")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <= PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45191")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <=> PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/20:0) + H+ + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:62537"))
    
    ## pe_to_nape_sn2
    pe <- "PE(14:0/16:0)"
    pc <- "PC(18:0/20:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45192")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) = PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 = CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45193")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) => PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 => CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45194")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <= PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <= CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:15378 + CHEBI:62537"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe, PC = pc), 
        template = list(), reaction = "RHEA:45195")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(18:0/20:0)")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(18:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PC(18:0/20:0) + PE(14:0/16:0) <=> PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PC(18:0/20:0) + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(18:0/0:0) + H+ + NAPE(14:0/16:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:57643 + CHEBI:64612 <=> CHEBI:58168 + CHEBI:15378 + CHEBI:62537")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:57643 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58168 + CHEBI:15378 + CHEBI:62537"))
    
    ## dhsm_to_dhcer
    dhsm <- "SM(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = list(), reaction =  "RHEA:45300")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O = Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 = CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64583 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = list(), reaction =  "RHEA:45301")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O => Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 => CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64583 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = list(), reaction =  "RHEA:45302")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O <= Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <= CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64583 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:15378 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(DhSM = dhsm), 
        template = list(), reaction =  "RHEA:45303")
    expect_equal(reaction_l[[1]]$DhSM, "SM(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O <=> Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("SM(d16:0(3OH,4OH)(15Me)/12:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + H+ + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64583 + CHEBI:15377 <=> CHEBI:83273 + CHEBI:15378 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64583 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:15378 + CHEBI:295975"))
    
    ## lnape_to_gpnae
    lnape <- "NAPE(14:0/0:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(LNAPE = lnape), 
        template = list(), reaction =  "RHEA:45420")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) = FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/0:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + GPNAE(0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 = CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:85216"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85225"))
    
    reaction_l <- .create_reaction(substrates = list(LNAPE = lnape), 
        template = list(), reaction =  "RHEA:45421")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) => FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/0:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + GPNAE(0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 => CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:85216"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85225"))
    
    reaction_l <- .create_reaction(substrates = list(LNAPE = lnape), 
        template = list(), reaction =  "RHEA:45422")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) <= FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/0:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + GPNAE(0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:85216"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85225"))
    
    reaction_l <- .create_reaction(substrates = list(LNAPE = lnape), 
        template = list(), reaction =  "RHEA:45423")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/0:0)")
    expect_equal(reaction_l[[1]]$GPNAE, "GPNAE(0:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/0:0/0:0) <=> FA(14:0) + H+ + GPNAE(0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/0:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + GPNAE(0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:85216 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85225")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:85216"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85225"))
    
    ## nape_to_lnape
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:45460")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) = FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(16:0) + H+ + NAPE(14:0/0:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 = CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85216"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:45461")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) => FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(16:0) + H+ + NAPE(14:0/0:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 => CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85216"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:45462")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <= FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(16:0) + H+ + NAPE(14:0/0:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <= CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85216"))
    
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "RHEA:45463")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$LNAPE, "NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + NAPE(14:0/16:0/18:0) <=> FA(16:0) + H+ + NAPE(14:0/0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + NAPE(14:0/16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(16:0) + H+ + NAPE(14:0/0:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:62537 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:85216")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:62537"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:85216"))
    
    ## sm_to_cer
    sm <- "SM(d16:1(3OH,4OH)(15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(SM = sm), 
        template = list(), reaction = "RHEA:45644")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(d16:1(3OH,4OH)(15Me)/12:0) = H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + SM(d16:1(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 = CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:78646"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15378 + CHEBI:72959 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(SM = sm), 
        template = list(), reaction = "RHEA:45645")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(d16:1(3OH,4OH)(15Me)/12:0) => H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + SM(d16:1(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 => CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:78646"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15378 + CHEBI:72959 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(SM = sm), 
        template = list(), reaction = "RHEA:45646")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(d16:1(3OH,4OH)(15Me)/12:0) <= H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + SM(d16:1(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <= CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:78646"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15378 + CHEBI:72959 + CHEBI:295975"))
    
    reaction_l <- .create_reaction(substrates = list(SM = sm), 
        template = list(), reaction = "RHEA:45647")
    expect_equal(reaction_l[[1]]$SM, "SM(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "H2O + SM(d16:1(3OH,4OH)(15Me)/12:0) <=> H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + SM(d16:1(3OH,4OH)(15Me)/12:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("H+ + Cer(d16:1(3OH,4OH)(15Me)/12:0) + Phosphocholine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:78646 <=> CHEBI:15378 + CHEBI:72959 + CHEBI:295975")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:78646"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:15378 + CHEBI:72959 + CHEBI:295975"))
    
    ## dhcer_to_cer 
    dhcer <- "Cer(d16:0(3OH,4OH)(15Me)/12:0)"
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:46544")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 = 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 = 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:46545")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 => 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 => 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:46546")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 <= 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <= 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377"))
    
    reaction_l <- .create_reaction(substrates = list(DhCer = dhcer),
        template = list(), reaction =  "RHEA:46547")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[1]]$Cer, "Cer(d16:1(3OH,4OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2 <=> 2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("Cer(d16:0(3OH,4OH)(15Me)/12:0) + 2 Fe2+-cytochrome b5 + 2 H+ + O2"))
    expect_equal(reaction_l[[2]]$reaction_product, c("2 Fe3+-cytochrome b5 + Cer(d16:1(3OH,4OH)(15Me)/12:0) + 2 H2O"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379 <=> 2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:31488 + 2 CHEBI:29033 + 2 CHEBI:15378 + CHEBI:15379"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("2 CHEBI:57540 + CHEBI:52639 + 2 CHEBI:15377"))
    
    ## coa_to_fao
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:52716")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH = FAO(18:0) + CoA + 2 NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + 2 H+ + 2 NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAO(18:0) + CoA + 2 NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 = CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:52717")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH => FAO(18:0) + CoA + 2 NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + 2 H+ + 2 NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAO(18:0) + CoA + 2 NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 => CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:52718")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH <= FAO(18:0) + CoA + 2 NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + 2 H+ + 2 NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAO(18:0) + CoA + 2 NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <= CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa),
        template = list(), reaction =  "RHEA:52719")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$FAO, "FAO(18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + 2 H+ + 2 NADPH <=> FAO(18:0) + CoA + 2 NADP+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + 2 H+ + 2 NADPH"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FAO(18:0) + CoA + 2 NADP+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783 <=> CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:83139 + 2 CHEBI:15378 + 2 CHEBI:57783"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:77396 + CHEBI:57287 + 2 CHEBI:58349"))
    
    ## sphinga_to_dhcer
    acylcoa <- "CoA(12:0)"
    sph <- "SPH(d16:0(1OH,3OH)(15Me))"
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa, SPH = sph), 
        template = list(), reaction = "RHEA:53424")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(1OH,3OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) = Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 = CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77636 + CHEBI:84410"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:57287 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa, SPH = sph), 
        template = list(), reaction = "RHEA:53425")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(1OH,3OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) => Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 => CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77636 + CHEBI:84410"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:57287 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa, SPH = sph), 
        template = list(), reaction = "RHEA:53426")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(1OH,3OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <= Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <= CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77636 + CHEBI:84410"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:57287 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(AcylCoA = acylcoa, SPH = sph), 
        template = list(), reaction = "RHEA:53427")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(12:0)")
    expect_equal(reaction_l[[1]]$DhCer, "Cer(d16:0(1OH,3OH)(15Me)/12:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me)) <=> Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(12:0) + SPH(d16:0(1OH,3OH)(15Me))"))
    expect_equal(reaction_l[[2]]$reaction_product, c("Cer(d16:0(1OH,3OH)(15Me)/12:0) + CoA + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77636 + CHEBI:84410 <=> CHEBI:83273 + CHEBI:57287 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77636 + CHEBI:84410"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:83273 + CHEBI:57287 + CHEBI:15378"))
    
    ## pep_to_napep_sn1
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = list(), reaction =  "RHEA:63596")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PE(P-16:0/14:0) + PC(20:0/18:0) = PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 = CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:140451"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = list(), reaction =  "RHEA:63597")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PE(P-16:0/14:0) + PC(20:0/18:0) => PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 => CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:140451"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = list(), reaction =  "RHEA:63598")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PE(P-16:0/14:0) + PC(20:0/18:0) <= PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <= CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:140451"))
    
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = list(), reaction =  "RHEA:63599")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PC(0:0/18:0) + H+ + NAPE(P-16:0/14:0/20:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:57875 + CHEBI:15378 + CHEBI:140451")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57875 + CHEBI:15378 + CHEBI:140451"))
    
    ## pe_to_dg
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:78951")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) = P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("P-Ethanolamine + DG(14:0/16:0/0:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 = CHEBI:58190 + CHEBI:17815 + CHEBI:15378") #########
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58190 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:78952")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) => P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("P-Ethanolamine + DG(14:0/16:0/0:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 => CHEBI:58190 + CHEBI:17815 + CHEBI:15378") #########
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58190 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:78953")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <= P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("P-Ethanolamine + DG(14:0/16:0/0:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <= CHEBI:58190 + CHEBI:17815 + CHEBI:15378") #########
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58190 + CHEBI:17815 + CHEBI:15378"))
    
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "RHEA:78954")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "H2O + PE(14:0/16:0) <=> P-Ethanolamine + DG(14:0/16:0/0:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("H2O + PE(14:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("P-Ethanolamine + DG(14:0/16:0/0:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:15377 + CHEBI:64612 <=> CHEBI:58190 + CHEBI:17815 + CHEBI:15378") #########
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:15377 + CHEBI:64612"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58190 + CHEBI:17815 + CHEBI:15378"))
    
    ## sn2lpe_to_fa
    sn2lpe <- "PE(0:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(sn2LPE = sn2lpe), 
        template = list(), reaction =  "sn2lpe_to_fa")
    expect_equal(reaction_l[[1]]$sn2LPE, "PE(0:0/14:0)")
    expect_equal(reaction_l[[1]]$FA, "FA(14:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(0:0/14:0) + H2O <=> FA(14:0) + H+ + Glycerophosphoethanolamine")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(0:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("FA(14:0) + H+ + Glycerophosphoethanolamine"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:65213 + CHEBI:15377 <=> CHEBI:28868 + CHEBI:15378 + CHEBI:143890")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:65213 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:28868 + CHEBI:15378 + CHEBI:143890"))
    
    ## lpeo_to_peo
    lpeo <- "PE(O-16:0/0:0)"
    acylcoa <- "CoA(18:0)"
    reaction_l <- .create_reaction(
        substrates = list(sn1LPEO = lpeo, AcylCoA = acylcoa), 
        template = list(), reaction =  "lpeo_to_peo")
    expect_equal(reaction_l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(reaction_l[[1]]$AcylCoA, "CoA(18:0)")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "CoA(18:0) + PE(O-16:0/0:0) <=> CoA + PE(O-16:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("CoA(18:0) + PE(O-16:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("CoA + PE(O-16:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:58342 +  <=> CHEBI:57287 + CHEBI:75028") ##############
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:58342 + "))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:57287 + CHEBI:75028"))
    
    ## sn2mg_to_sn1mg
    sn2mg <- "MG(0:0/14:0/0:0)"
    reaction_l <- .create_reaction(substrates = list(sn2MG = sn2mg), 
        template = list(), reaction =  "sn2mg_to_sn1mg")
    expect_equal(reaction_l[[1]]$sn2MG, "MG(0:0/14:0/0:0)")
    expect_equal(reaction_l[[1]]$sn1MG, "MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "MG(0:0/14:0/0:0) <=> MG(14:0/0:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("MG(0:0/14:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("MG(14:0/0:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:35759 <=> CHEBI:17389")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:35759"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:17389"))
    
    ## nape_to_pnae
    nape <- "NAPE(14:0/16:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(NAPE = nape), 
        template = list(), reaction =  "nape_to_pnae")
    expect_equal(reaction_l[[1]]$NAPE, "NAPE(14:0/16:0/18:0)")
    expect_equal(reaction_l[[1]]$PNAE, "PNAE(18:0)")
    expect_equal(reaction_l[[1]]$DG, "DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(14:0/16:0/18:0) + H2O <=> PNAE(18:0) + DG(14:0/16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("NAPE(14:0/16:0/18:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PNAE(18:0) + DG(14:0/16:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:62537 + CHEBI:15377 <=> CHEBI:145538 + CHEBI:17815")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:62537 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:145538 + CHEBI:17815"))
    
    ## napeo_to_nae
    napeo <- "NAPE(O-18:0/16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(NAPEO = napeo), 
        template = list(), reaction =  "napeo_to_nae")
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-18:0/16:0/14:0)")
    expect_equal(reaction_l[[1]]$NAE, "NAE(14:0)")
    expect_equal(reaction_l[[1]]$PAO, "PA(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "NAPE(O-18:0/16:0/14:0) + H2O <=> NAE(14:0) + PA(O-18:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("NAPE(O-18:0/16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("NAE(14:0) + PA(O-18:0/16:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, " + CHEBI:15377 <=> CHEBI:52640 + CHEBI:73332") ###########
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c(" + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:52640 + CHEBI:73332"))

    ## pe_to_pa
    pe <- "PE(14:0/16:0)"
    reaction_l <- .create_reaction(substrates = list(PE = pe), 
        template = list(), reaction = "pe_to_pa")
    expect_equal(reaction_l[[1]]$PE, "PE(14:0/16:0)")
    expect_equal(reaction_l[[1]]$PA, "PA(14:0/16:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(14:0/16:0) + H2O <=> PA(14:0/16:0) + Ethanolamine + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(14:0/16:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PA(14:0/16:0) + Ethanolamine + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:64612 + CHEBI:15377 <=> CHEBI:58608 + CHEBI:57603 + CHEBI:15378")
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:64612 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:58608 + CHEBI:57603 + CHEBI:15378"))
    
    ## peo_to_lpeo
    peo <- "PE(O-16:0/14:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo), 
        template = list(), reaction =  "peo_to_lpeo")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$sn1LPEO, "PE(O-16:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + H2O <=> PE(O-16:0/0:0) + FA(14:0) + H+")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + H2O"))
    expect_equal(reaction_l[[2]]$reaction_product, c("PE(O-16:0/0:0) + FA(14:0) + H+"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:15377 <=>  + CHEBI:28868 + CHEBI:15378") ######
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:15377"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c(" + CHEBI:28868 + CHEBI:15378"))
    
    ## peo_to_napeo_sn1
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = list(), reaction =  "peo_to_napeo_sn1") ################
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/20:0)")
    expect_equal(reaction_l[[1]]$sn2LPC, "PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("NAPE(O-16:0/14:0/20:0) + PC(0:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ######
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c(" + CHEBI:58168"))
    
    ## peo_to_napeo_sn2
    peo <- "PE(O-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEO = peo, PC = pc), 
        template = list(), reaction =  "peo_to_napeo_sn2")
    expect_equal(reaction_l[[1]]$PEO, "PE(O-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEO, "NAPE(O-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_name, "")
    expect_equal(reaction_l[[2]]$reaction_formula, 
        "PE(O-16:0/14:0) + PC(20:0/18:0) <=> NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(O-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("NAPE(O-16:0/14:0/18:0) + PC(20:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:75028 + CHEBI:57643 <=>  + CHEBI:58168") ######
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:75028 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c(" + CHEBI:58168"))
    
    ## pep_to_napep_sn2
    pep <- "PE(P-16:0/14:0)"
    pc <- "PC(20:0/18:0)"
    reaction_l <- .create_reaction(substrates = list(PEP = pep, PC = pc), 
        template = list(), reaction =  "pep_to_napep_sn2")
    expect_equal(reaction_l[[1]]$PEP, "PE(P-16:0/14:0)")
    expect_equal(reaction_l[[1]]$PC, "PC(20:0/18:0)")
    expect_equal(reaction_l[[1]]$NAPEP, "NAPE(P-16:0/14:0/18:0)")
    expect_equal(reaction_l[[1]]$sn1LPC, "PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_formula,
        "PE(P-16:0/14:0) + PC(20:0/18:0) <=> NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)")
    expect_equal(reaction_l[[2]]$reaction_isReversible, "")
    expect_equal(reaction_l[[2]]$reaction_geneAssociation, "")
    expect_equal(reaction_l[[2]]$reaction_pathway, "")
    expect_equal(reaction_l[[2]]$reaction_substrate, c("PE(P-16:0/14:0) + PC(20:0/18:0)"))
    expect_equal(reaction_l[[2]]$reaction_product, c("NAPE(P-16:0/14:0/18:0) + PC(20:0/0:0)"))
    expect_equal(reaction_l[[2]]$reaction_formula_chebi, "CHEBI:77290 + CHEBI:57643 <=> CHEBI:140451 + CHEBI:58168") ######
    expect_equal(reaction_l[[2]]$reaction_substrate_chebi, c("CHEBI:77290 + CHEBI:57643"))
    expect_equal(reaction_l[[2]]$reaction_product_chebi, c("CHEBI:140451 + CHEBI:58168"))
    
})


## function create_reactions
test_that("create_reactions works.", {
    FA <- c("FA(16:0)", "FA(12:0)", "FA(14:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "", "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE, "unknown", "gene_1", "pathway_1"),
        c(2, "", "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE, "unknown", "gene_2", "pathway_2"),
        c(3, "", "RHEA:19709", "M_LPA + M_AcylCoA = M_CoA + M_PA", FALSE, "unknown", "gene_3", "pathway_3"),
        c(4, "", "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE, "unknown", "gene_4", "pathway_4"))
    
    reactions <- data.frame(order = reactions[, 1], 
        reaction_name = reactions[, 2], reaction_RHEA = reactions[, 3],
        reaction_formula = reactions[, 4], directed = reactions[, 5],
        reaction_isReversible = reactions[, 6],
        reaction_geneAssociation = reactions[, 7],
        reaction_pathway = reactions[, 8])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
    
    ## run the function
    reactions_l <- create_reactions(substrates = list(FA = FA), 
        reactions = reactions)
    expect_equal(length(reactions_l), 4)
    
    ## first entry
    expect_equal(names(reactions_l[[1]][[1]]), c("FA", "AcylCoA"))
    expect_equal(reactions_l[[1]][[1]]$FA, 
        c("FA(16:0)", "FA(12:0)", "FA(14:0)"))
    expect_equal(reactions_l[[1]][[1]]$AcylCoA, 
        c("CoA(16:0)", "CoA(12:0)", "CoA(14:0)"))
    expect_equal(names(reactions_l[[1]][[2]]), 
        c("order", "reaction_name", "reaction_RHEA", "reaction_formula", 
            "directed", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_constraints", "reaction_constraints_negate", 
            "reaction_substrate", "reaction_product", 
            "reaction_formula_chebi", "reaction_substrate_chebi", 
            "reaction_product_chebi"))
    expect_equal(length(reactions_l[[1]][[2]]$reaction_formula), 3)
    expect_equal(reactions_l[[1]][[2]]$reaction_formula,
        c("ATP + CoA + FA(16:0) = PPi + AMP + CoA(16:0)",
            "ATP + CoA + FA(12:0) = PPi + AMP + CoA(12:0)",
            "ATP + CoA + FA(14:0) = PPi + AMP + CoA(14:0)"))
    expect_equal(reactions_l[[1]][[2]]$reaction_isReversible, 
        "unknown")
    expect_equal(reactions_l[[1]][[2]]$reaction_geneAssociation,
        "gene_1")
    expect_equal(reactions_l[[1]][[2]]$reaction_pathway,
        "pathway_1")
    expect_equal(reactions_l[[1]][[2]]$reaction_substrate,
        c("ATP + CoA + FA(16:0)", "ATP + CoA + FA(12:0)", "ATP + CoA + FA(14:0)"))
    expect_equal(reactions_l[[1]][[2]]$reaction_product,
        c("PPi + AMP + CoA(16:0)", "PPi + AMP + CoA(12:0)", "PPi + AMP + CoA(14:0)"))
    
    ## second entry
    expect_equal(names(reactions_l[[2]][[1]]), c("AcylCoA", "sn1LPA"))
    expect_equal(reactions_l[[2]][[1]]$AcylCoA, 
        c("CoA(16:0)", "CoA(12:0)", "CoA(14:0)"))
    expect_equal(reactions_l[[2]][[1]]$sn1LPA, 
        c("PA(16:0/0:0)", "PA(12:0/0:0)", "PA(14:0/0:0)"))
    expect_equal(names(reactions_l[[2]][[2]]), 
        c("order", "reaction_name", "reaction_RHEA", "reaction_formula", 
            "directed", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_constraints", "reaction_constraints_negate", 
            "reaction_substrate", "reaction_product", 
            "reaction_formula_chebi", "reaction_substrate_chebi", 
            "reaction_product_chebi"))
    expect_equal(length(reactions_l[[2]][[2]]$reaction_formula), 3)
    expect_equal(reactions_l[[2]][[2]]$reaction_formula,
        c("Glycerol-3-P + CoA(16:0) = CoA + PA(16:0/0:0)",
            "Glycerol-3-P + CoA(12:0) = CoA + PA(12:0/0:0)",
            "Glycerol-3-P + CoA(14:0) = CoA + PA(14:0/0:0)"))
    expect_equal(reactions_l[[2]][[2]]$reaction_isReversible, 
        "unknown")
    expect_equal(reactions_l[[2]][[2]]$reaction_geneAssociation, 
        "gene_2")
    expect_equal(reactions_l[[2]][[2]]$reaction_pathway, 
        "pathway_2")
    expect_equal(reactions_l[[2]][[2]]$reaction_substrate,
        c("Glycerol-3-P + CoA(16:0)", "Glycerol-3-P + CoA(12:0)",
            "Glycerol-3-P + CoA(14:0)"))
    expect_equal(reactions_l[[2]][[2]]$reaction_product,
        c("CoA + PA(16:0/0:0)", "CoA + PA(12:0/0:0)", "CoA + PA(14:0/0:0)"))
    
    ## third entry
    expect_equal(names(reactions_l[[3]][[1]]), c("sn1LPA", "AcylCoA", "PA"))
    expect_equal(reactions_l[[3]][[1]]$sn1LPA, 
        c("PA(16:0/0:0)", "PA(12:0/0:0)", "PA(14:0/0:0)"))
    expect_equal(reactions_l[[3]][[1]]$AcylCoA, 
        c("CoA(16:0)", "CoA(12:0)", "CoA(14:0)"))
    expect_equal(reactions_l[[3]][[1]]$PA, 
        c("PA(16:0/16:0)", "PA(12:0/16:0)", "PA(14:0/16:0)", "PA(16:0/12:0)",
            "PA(12:0/12:0)", "PA(14:0/12:0)", "PA(16:0/14:0)", "PA(12:0/14:0)",
            "PA(14:0/14:0)"))
    expect_equal(names(reactions_l[[3]][[2]]), 
          c("order", "reaction_name", "reaction_RHEA", "reaction_formula", 
            "directed", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_constraints", "reaction_constraints_negate", 
            "reaction_substrate", "reaction_product", 
            "reaction_formula_chebi", "reaction_substrate_chebi", 
            "reaction_product_chebi"))
    expect_equal(length(reactions_l[[3]][[2]]$reaction_formula), 9)
    expect_equal(reactions_l[[3]][[2]]$reaction_formula,
        c("PA(16:0/0:0) + CoA(16:0) = CoA + PA(16:0/16:0)",
            "PA(12:0/0:0) + CoA(16:0) = CoA + PA(12:0/16:0)",
            "PA(14:0/0:0) + CoA(16:0) = CoA + PA(14:0/16:0)", 
            "PA(16:0/0:0) + CoA(12:0) = CoA + PA(16:0/12:0)",
            "PA(12:0/0:0) + CoA(12:0) = CoA + PA(12:0/12:0)", 
            "PA(14:0/0:0) + CoA(12:0) = CoA + PA(14:0/12:0)",
            "PA(16:0/0:0) + CoA(14:0) = CoA + PA(16:0/14:0)", 
            "PA(12:0/0:0) + CoA(14:0) = CoA + PA(12:0/14:0)",
            "PA(14:0/0:0) + CoA(14:0) = CoA + PA(14:0/14:0)"))
    expect_equal(reactions_l[[3]][[2]]$reaction_isReversible, 
        "unknown")
    expect_equal(reactions_l[[3]][[2]]$reaction_geneAssociation, 
        "gene_3")
    expect_equal(reactions_l[[3]][[2]]$reaction_pathway, 
        "pathway_3")
    expect_equal(reactions_l[[3]][[2]]$reaction_substrate,
        c("PA(16:0/0:0) + CoA(16:0)", "PA(12:0/0:0) + CoA(16:0)",
            "PA(14:0/0:0) + CoA(16:0)", "PA(16:0/0:0) + CoA(12:0)",
            "PA(12:0/0:0) + CoA(12:0)", "PA(14:0/0:0) + CoA(12:0)",
            "PA(16:0/0:0) + CoA(14:0)", "PA(12:0/0:0) + CoA(14:0)",
            "PA(14:0/0:0) + CoA(14:0)"))
    expect_equal(reactions_l[[3]][[2]]$reaction_product,
        c("CoA + PA(16:0/16:0)", "CoA + PA(12:0/16:0)", "CoA + PA(14:0/16:0)",
            "CoA + PA(16:0/12:0)", "CoA + PA(12:0/12:0)", "CoA + PA(14:0/12:0)",
            "CoA + PA(16:0/14:0)", "CoA + PA(12:0/14:0)", "CoA + PA(14:0/14:0)"))
    
    ## fourth entry
    expect_equal(names(reactions_l[[4]][[1]]), c("PA", "DG"))
    expect_equal(reactions_l[[4]][[1]]$PA, 
        c("PA(16:0/16:0)", "PA(12:0/16:0)", "PA(14:0/16:0)", "PA(16:0/12:0)",
            "PA(12:0/12:0)", "PA(14:0/12:0)", "PA(16:0/14:0)", "PA(12:0/14:0)",
            "PA(14:0/14:0)"))
    expect_equal(reactions_l[[4]][[1]]$DG, 
        c("DG(16:0/16:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/16:0/0:0)",
            "DG(16:0/12:0/0:0)", "DG(12:0/12:0/0:0)", "DG(14:0/12:0/0:0)",
            "DG(16:0/14:0/0:0)", "DG(12:0/14:0/0:0)", "DG(14:0/14:0/0:0)"))
    expect_equal(names(reactions_l[[4]][[2]]), 
       c("order", "reaction_name", "reaction_RHEA", "reaction_formula", 
            "directed", "reaction_isReversible",
            "reaction_geneAssociation", "reaction_pathway", 
            "reaction_constraints", "reaction_constraints_negate", 
            "reaction_substrate", "reaction_product", 
            "reaction_formula_chebi", "reaction_substrate_chebi", 
            "reaction_product_chebi"))
    expect_equal(length(reactions_l[[1]][[2]]$reaction_formula), 3)
    expect_equal(reactions_l[[4]][[2]]$reaction_formula,
        c("H2O + PA(16:0/16:0) = Pi + DG(16:0/16:0/0:0)",
          "H2O + PA(12:0/16:0) = Pi + DG(12:0/16:0/0:0)",
          "H2O + PA(14:0/16:0) = Pi + DG(14:0/16:0/0:0)",
          "H2O + PA(16:0/12:0) = Pi + DG(16:0/12:0/0:0)",
          "H2O + PA(12:0/12:0) = Pi + DG(12:0/12:0/0:0)",
          "H2O + PA(14:0/12:0) = Pi + DG(14:0/12:0/0:0)",
          "H2O + PA(16:0/14:0) = Pi + DG(16:0/14:0/0:0)",
          "H2O + PA(12:0/14:0) = Pi + DG(12:0/14:0/0:0)",
          "H2O + PA(14:0/14:0) = Pi + DG(14:0/14:0/0:0)"))
    expect_equal(reactions_l[[4]][[2]]$reaction_isReversible, 
        "unknown")
    expect_equal(reactions_l[[4]][[2]]$reaction_geneAssociation, 
        "gene_4")
    expect_equal(reactions_l[[4]][[2]]$reaction_pathway, 
        "pathway_4")
    expect_equal(reactions_l[[4]][[2]]$reaction_substrate,
        c("H2O + PA(16:0/16:0)", "H2O + PA(12:0/16:0)", "H2O + PA(14:0/16:0)",
            "H2O + PA(16:0/12:0)", "H2O + PA(12:0/12:0)", "H2O + PA(14:0/12:0)",
            "H2O + PA(16:0/14:0)", "H2O + PA(12:0/14:0)", "H2O + PA(14:0/14:0)"))
    expect_equal(reactions_l[[4]][[2]]$reaction_product,
        c("Pi + DG(16:0/16:0/0:0)", "Pi + DG(12:0/16:0/0:0)",
          "Pi + DG(14:0/16:0/0:0)", "Pi + DG(16:0/12:0/0:0)",
          "Pi + DG(12:0/12:0/0:0)", "Pi + DG(14:0/12:0/0:0)",
          "Pi + DG(16:0/14:0/0:0)", "Pi + DG(12:0/14:0/0:0)",
          "Pi + DG(14:0/14:0/0:0)"))
})

## function create_reaction_adjacency_matrix
test_that("create_reaction_adjacency_matrix works.", {
    FA <- c("FA(16:0)", "FA(12:0)", "FA(14:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "", "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE, "unknown", "gene_1", "pathway_1"),
        c(2, "", "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE, "unknown", "gene_2", "pathway_2"),
        c(3, "", "RHEA:19709", "M_LPA + M_AcylCoA = M_CoA + M_PA", FALSE, "unknown", "gene_3", "pathway_3"),
        c(4, "", "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE, "unknown", "gene_4", "pathway_4"),
        c(5, "", "RHEA:52716", "M_AcylCoA + 2 M_NADPH + 2 M_H+ = M_FAO + 2 M_NADP + M_CoA", FALSE, "unknown", "gene_5", "pathway_5"))
    
    reactions <- data.frame(order = reactions[, 1], 
        reaction_name = reactions[, 2], reaction_RHEA = reactions[, 3],
        reaction_formula = reactions[, 4], directed = reactions[, 5],
        reaction_isReversible = reactions[, 6],
        reaction_geneAssociation = reactions[, 7],
        reaction_pathway = reactions[, 8])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
    
    reactions_l <- create_reactions(substrates = list(FA = FA), 
        reactions = reactions)
    
    ## run the function
    adj <- create_reaction_adjacency_matrix(reaction_l = reactions_l)
    expect_equal(dim(adj), c(40, 40))
    expect_equal(rownames(adj), colnames(adj))
    expect_equal(sort(rownames(adj))[1:38], ## missed Pi and PPi
        c("AMP", "ATP", "CoA", "CoA(12:0)", "CoA(14:0)", "CoA(16:0)", 
            "DG(12:0/12:0/0:0)", "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)",
            "DG(14:0/12:0/0:0)", "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)",
            "DG(16:0/12:0/0:0)", "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)",
            "FA(12:0)", "FA(14:0)", "FA(16:0)", "FAO(12:0)", "FAO(14:0)", 
            "FAO(16:0)", "Glycerol-3-P", "H+", "H2O", "NADP+", "NADPH",
            "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)", "PA(12:0/16:0)",
            "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", "PA(14:0/16:0)",
            "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)", "PA(16:0/16:0)",
            "Pi", "PPi")[1:38])
    expect_equal(sum(adj), 94)
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


