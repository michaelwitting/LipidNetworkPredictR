## function .create_reaction
test_that(".check_reaction works", {
    
    ## throws error as expected
    expect_error(LipidNetworkPredictR:::.check_reaction(
        reaction = c("RHEA:36171", "RHEA:36172")),
        "'reaction' has to be of length 1")
    
    expect_error(LipidNetworkPredictR:::.check_reaction(reaction = "foo"),
        "foo is not an implemented reaction type")
    expect_error(LipidNetworkPredictR:::.check_reaction(reaction = "RHEA:"),
        "RHEA: is not an implemented reaction type")
    expect_error(LipidNetworkPredictR:::.check_reaction(reaction = 1),
        "1 is not an implemented reaction type")
    
    
})
