## function .create_reaction
test_that(".create_substrates_combination works", {
    
    acyldhap <- c("DHAP(18:0)", "DHAP(20:0)")
    fao <- c("FAO(16:0)", "FAO(18:0)", "FAO(20:0)") 
    substrates <- list(AcylDHAP = acyldhap, FAO = fao)
    
    reaction <- "RHEA:36171"
    template <- .create_template(template = list(), reaction = reaction)
    df_substrates <- .create_substrates_combinations(
        substrates = substrates, template = template)
    expect_equal(df_substrates$AcylDHAP, rep(c("DHAP(18:0)", "DHAP(20:0)"), 3))
    expect_equal(df_substrates$FAO, c("FAO(16:0)", "FAO(16:0)", 
        "FAO(18:0)", "FAO(18:0)", "FAO(20:0)", "FAO(20:0)"))
    
    substrates <- list(AcylDHAP = NULL, FAO = NULL)
    expect_error(.create_substrates_combinations(
        substrates = substrates, template = template), 
        "must be the same length as the vector")
   
})
