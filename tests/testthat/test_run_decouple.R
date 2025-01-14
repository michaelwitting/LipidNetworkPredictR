## function run_redouple
test_that("run_decouple works", {
    
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
     
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
     
    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
    
    ## create scores (simulate scores via rnorm with sd = 2)
    scores <- rbind(
        c("AMP", 0.94446911),
        c("ATP", 0.78505230),
        c("CoA", 2.90333236),
        c("CoA(12:0)", 0.13577818),
        c("CoA(14:0)", 0.64465176),
        c("CoA(16:0)", 0.26128554),
        c("FA(12:0)", -1.98870069),
        c("FA(14:0)", -2.17736671),
        c("FA(16:0)", -2.41113650),
        c("Glycerol-3-P", -0.84185713 ),
        c("PA(12:0/0:0)", -2.12108629),
        c("PA(12:0/12:0)", 1.60395612),
        c("PA(12:0/14:0)", 0.07693279),
        c("PA(12:0/16:0)", 2.46639841),
        c("PA(14:0/0:0)", -1.73630919),
        c("PA(14:0/12:0)", 1.32713790),
        c("PA(14:0/14:0)", 3.78234917),
        c("PA(14:0/16:0)", 1.88168721),
        c("PA(16:0/0:0)", -0.84712629),
        c("PA(16:0/12:0)", 4.16900272),
        c("PA(16:0/14:0)", 2.07352961),
        c("PA(16:0/16:0)", 4.48556057),
        c("PPi", -0.28139118)
    )
    scores <- data.frame(feature = scores[, 1], score = scores[, 2])
    scores$feature <- as.character(scores$feature)
    scores$score <- as.numeric(scores$score)

    ## run the function
    scores_default <- run_decouple(reaction_l = reaction_l, scores = scores)
    
    expect_equal(dim(scores_default), c(32, 6))
    expect_equal(scores_default[["run_id"]], 
        c(rep(1, 4), rep(2, 4), rep(3, 12), rep(4, 12)))
    expect_equal(scores_default[["statistic"]], 
        rep(c("mlm", "ulm", "corr_wsum", "norm_wsum", "wsum", "corr_wmean",
            "norm_wmean", "wmean"), each = 4))
    expect_equal(unique(scores_default[["source"]]), 
        c("CHEBI:58342", "CHEBI:57560", "CHEBI:57970", "CHEBI:58608"))
    expect_equal(scores_default[["condition"]], rep("score", 32))
    expect_equal(scores_default[["score"]], 
        c(-0.4018735, -3.2794148, -2.5721334,  2.5630578, -2.9855025, 
            -2.1597864, -0.2732775,  4.5209298, -4.8982246, -1.7724319,
            0.4145403, 37.1506202, -2.4035082, -2.2962217, -0.2798554,
            3.3275044, -6.5772039, -4.7045218, 1.0417155, 21.8665545,
            -1.6327415, -0.5908106, 0.1381801, 4.1278467, -2.4035082,
            -2.2962217, -0.2798554, 3.3275044, -2.1924013, -1.5681739,
            0.3472385, 2.4296172),
        tolerance = 1e-06)
    expect_equal(scores_default[["p_value"]], 
        c(0.6925087653, 0.0041663126, 0.0191897122, 0.0195580501, 0.0070514841,
            0.0425031401, 0.7873094619, 0.0001870385, 0.18, 0.42, 0.4, 0.02,
            0.18, 0.42, 0.4, 0.02, 0.18, 0.42, 0.4, 0.02, 0.18, 0.42, 0.4,
            0.02, 0.18, 0.42, 0.4, 0.02, 0.18, 0.42, 0.4, 0.02),
        tolerance = 1e-06)
    
    ## times = 10
    scores_times10 <- run_decouple(reaction_l = reaction_l, scores = scores, times = 10)
    expect_equal(dim(scores_times10), c(32, 6))
    expect_equal(scores_times10[["run_id"]], 
        c(rep(1, 4), rep(2, 4), rep(3, 12), rep(4, 12)))
    expect_equal(scores_times10[["statistic"]], 
        rep(c("mlm", "ulm", "corr_wsum", "norm_wsum", "wsum", "corr_wmean",
            "norm_wmean", "wmean"), each = 4))
    expect_equal(unique(scores_times10[["source"]]), 
        c("CHEBI:58342", "CHEBI:57560", "CHEBI:57970", "CHEBI:58608"))
    expect_equal(scores_times10[["condition"]], rep("score", 32))
    expect_equal(scores_times10[["score"]], 
        c(-0.40187347, -3.27941484, -2.57213340, 2.56305782, -2.98550253,
            -2.15978641, -0.27327751, 4.52092985, -4.59726824, -3.28831960, 0,
            15.28406569, -2.24118890, -2.80769378, -0.04151253, 3.03051578,
            -6.57720390, -4.70452177, 1.04171548, 21.86655450, -1.53242275,
            -1.09610653, 0, 1.69822952, -2.24118890, -2.80769378, -0.04151253,
            3.03051578, -2.19240130, -1.56817392, 0.34723849, 2.42961717),
        tolerance = 1e-06)
    expect_equal(scores_times10[["p_value"]], 
        c(0.6925087653, 0.0041663126, 0.0191897122, 0.0195580501, 0.0070514841,
            0.0425031401, 0.7873094619, 0.0001870385, 0.2, 0.2, 1, 0.2, 0.2,
            0.2, 1, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1,
            0.2, 0.2, 0.2, 1, 0.2),
        tolerance = 1e-06)
    
    ## truncated reaction_l
    scores_truncated <- run_decouple(reaction_l = reaction_l[1:2], 
        scores = scores)
    expect_equal(dim(scores_truncated), c(24, 6))
    expect_equal(scores_truncated[["run_id"]], 
        c(rep(1, 3), rep(2, 3), rep(3, 9), rep(4, 9)))
    expect_equal(scores_truncated[["statistic"]], 
        rep(c("mlm", "ulm", "corr_wsum", "norm_wsum", "wsum", "corr_wmean",
            "norm_wmean", "wmean"), each = 3))
    expect_equal(unique(scores_truncated[["source"]]), 
                 c("CHEBI:58342", "CHEBI:57560", "CHEBI:57970"))
    expect_equal(scores_truncated[["condition"]], rep("score", 24))
    expect_equal(scores_truncated[["score"]], 
        c(-1.6759822, -4.5806860, -3.8667282, -2.9855025, -2.1597864,
            -0.2732775, -4.8982246, -1.7724319, 0.4145403, -2.4035082,
            -2.2962217, -0.2798554, -6.5772039, -4.7045218, 1.0417155,
            -1.6327415, -0.5908106, 0.1381801, -2.4035082, -2.2962217,
            -0.2798554, -2.1924013, -1.5681739, 0.3472385),
        tolerance = 1e-06)
    expect_equal(scores_truncated[["p_value"]], 
        c(0.1101199335, 0.0002042001, 0.0010388217, 0.0070514841, 0.0425031401,
            0.7873094619, 0.18, 0.42, 0.4, 0.18, 0.42, 0.4, 0.18, 0.42, 0.4, 
            0.18, 0.42, 0.4, 0.18, 0.42, 0.4, 0.18, 0.42, 0.4),
        tolerance = 1e-06)
    
    ## truncated scores
    scores_truncated <- run_decouple(reaction_l = reaction_l, 
        scores = scores[1:10, ])
    expect_equal(dim(scores_truncated), c(16, 6))
    expect_equal(scores_truncated[["run_id"]], 
        c(rep(1, 2), rep(2, 2), rep(3, 6), rep(4, 6)))
    expect_equal(scores_truncated[["statistic"]], 
        rep(c("mlm", "ulm", "corr_wsum", "norm_wsum", "wsum", "corr_wmean",
            "norm_wmean", "wmean"), each = 2))
    expect_equal(unique(scores_truncated[["source"]]), 
        c("CHEBI:58342", "CHEBI:57560"))
    expect_equal(scores_truncated[["condition"]], rep("score", 16))
    expect_equal(scores_truncated[["score"]], 
        c(-0.77079394, -4.03058444, -4.20343755, 0.62111591, -11.17447214,
            0.27876956, -2.38077150, 0.68057973, -6.57720390, 1.04171548, 
            -3.72482405, 0.09292319, -2.38077150, 0.68057973, -2.19240130,
            0.34723849),
        tolerance = 1e-06)
    expect_equal(scores_truncated[["p_value"]], 
        c(0.466044085, 0.004992096, 0.002982475, 0.551805525, 0.02, 0.54, 0.02,
            0.54, 0.02, 0.54, 0.02, 0.54, 0.02, 0.54, 0.02, 0.54),
        tolerance = 1e-06)
})

## function plot_scores
test_that("plot_scores works", {
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")

    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)

    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)

    ## create scores (simulate scores via rnorm with sd = 2)
    scores <- rbind(
        c("AMP", 0.94446911),
        c("ATP", 0.78505230),
        c("CoA", 2.90333236),
        c("CoA(12:0)", 0.13577818),
        c("CoA(14:0)", 0.64465176),
        c("CoA(16:0)", 0.26128554),
        c("FA(12:0)", -1.98870069),
        c("FA(14:0)", -2.17736671),
        c("FA(16:0)", -2.41113650),
        c("Glycerol-3-P", -0.84185713 ),
        c("PA(12:0/0:0)", -2.12108629),
        c("PA(12:0/12:0)", 1.60395612),
        c("PA(12:0/14:0)", 0.07693279),
        c("PA(12:0/16:0)", 2.46639841),
        c("PA(14:0/0:0)", -1.73630919),
        c("PA(14:0/12:0)", 1.32713790),
        c("PA(14:0/14:0)", 3.78234917),
        c("PA(14:0/16:0)", 1.88168721),
        c("PA(16:0/0:0)", -0.84712629),
        c("PA(16:0/12:0)", 4.16900272),
        c("PA(16:0/14:0)", 2.07352961),
        c("PA(16:0/16:0)", 4.48556057),
        c("PPi", -0.28139118)
    )
    scores <- data.frame(feature = scores[, 1], score = scores[, 2])
    scores$feature <- as.character(scores$feature)
    scores$score <- as.numeric(scores$score)
    
    ## run run_decouple
    scores <- run_decouple(reaction_l = reaction_l, scores = scores)
    
    ## run the function
    g <- plot_scores(scores)
    expect_equal(g, "gg")
    expect_error(plot_scores(scores = NULL), "column 'run_id' not in 'scores'")
    expect_error(plot_scores(scores = scores[, 1, drop = FALSE]),
        "column 'statistic' not in 'scores'")
    expect_error(plot_scores(scores = scores[, 1:2, drop = FALSE]),
        "column 'source' not in 'scores'")
    expect_error(plot_scores(scores = scores[, 1:3, drop = FALSE]),
        "column 'condition' not in 'scores'")
    expect_error(plot_scores(scores = scores[, 1:4, drop = FALSE]),
        "column 'object' not in 'scores'")
    expect_error(plot_scores(scores = scores[, 1:5, drop = FALSE]),
        "column 'p_value' not in 'scores'")
})
