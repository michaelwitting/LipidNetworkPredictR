## function add_attributes
test_that("add_attributes works", {
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
     
    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
     
    ## create the adjacency matrix
    reaction_adj <- create_reaction_adjacency_matrix(reaction_l)
    
    g <- igraph::graph_from_adjacency_matrix(reaction_adj, weighted = TRUE, 
        diag = FALSE)
    
    expect_error(
        add_attributes(g, attribute_type = "foo", attributes = NULL, cols_vertex = "name"), 
        "'arg' should be one of ")

    ## attribute_type: vertex
    attributes_df <- data.frame(
        name = c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)", "DG(12:0/12:0/0:0)",
            "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/12:0/0:0)",
            "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)", "DG(16:0/12:0/0:0)",
            "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)", "FA(12:0)", "FA(14:0)",
            "FA(16:0)", "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)",
            "PA(12:0/16:0)", "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", 
            "PA(14:0/16:0)", "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)",
            "PA(16:0/16:0)"),
        logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
            -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
            -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
            -0.11),
        logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
        -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
        -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
  
    ## apply the function
    g_new <- add_attributes(g, attribute_type = "vertex", attributes = attributes_df, 
        cols_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        c(NA, NA, NA, -5.08, 0.75, 5.43, -0.62, 2.35, 1.39, 2.91, 0.26, -4.14,
            0.19, 6.18, 0.78, -1.81, 4.66, -0.10, NA, NA, 2.84, -0.81, -0.81,
            -0.32, 0.17, 2.25, -1.94, 0.80, 4.21, 0.20, -3.29, -0.11, NA, NA))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        c(NA, NA, NA, -2.73, 6.14, 1.98, 0.09, 1.57, 1.77, 3.08, 4.04, -3.01,
            1.22, -4.25, 0.39, 0.53, 3.30, 7.10, NA, NA, 2.81, -0.99, -0.09,
            -8.25, 4.94, -3.54, -7.74, -1.98, 0.73, 2.36, 2.53, -0.62, NA, NA))
    expect_error(
        add_attributes(g, attribute_type = "vertex", attributes = attributes_df, 
            cols_vertex = "foo"), 
        "not in ")

    ## attribute_type: edges, attributes: data.frame
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
            c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
    ))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## some vertex names not in g
    attributes_df <- data.frame(
        name = c("foo_1", "foo_2", "CoA(16:0)"),
        logFC_cond1 = c(-5.08,  0.75,  5.43),
        logFC_cond2 = c(-2.73,  6.14,  1.98))
    g_new <- add_attributes(g, attribute_type = "vertex", 
        attributes = attributes_df, cols_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        c(NA, NA, NA, NA, NA, 5.43, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        c(NA, NA, NA, NA, NA, 1.98, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    
    ## all vertex names not in g
    attributes_df <- data.frame(
        name = c("foo_1", "foo_2", "foo_3"),
        logFC_cond1 = c(-5.08,  0.75,  5.43),
        logFC_cond2 = c(-2.73,  6.14,  1.98))
    g_new <- add_attributes(g, attribute_type = "vertex", 
        attributes = attributes_df, cols_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA)))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA)))
    
    ## apply the function
    g_new <- add_attributes(g, attribute_type = "edges", 
        attributes = attributes, cols_vertex = c("from", "to"))
    expect_equal(igraph::E(g_new)$weight, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, 0.8, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$weight[12], 0.5)
    expect_equal(attr(igraph::E(g_new), "vnames")[14], "CoA(12:0)|PA(14:0/12:0)")
    expect_equal(igraph::E(g_new)$weight[14], 0.8)
    expect_error(
        add_attributes(g, attribute_type = "edges", attributes = attributes, 
            cols_vertex = "name"), 
        " not in ")
    expect_error(
        add_attributes(g, attribute_type = "edges", attributes = attributes, 
            cols_vertex = c("from", "to", "weight")), 
        "has to be of length 2")

    ## attribute_type: edges, attributes: matrix
    attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
        dimnames = list(
        c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)"),
        c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)")))
    
    ## apply the function
    g_new <- add_attributes(g, attribute_type = "edges", attributes = attributes)
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$value[12], 0.5)
    expect_equal(attr(igraph::E(g_new), "vnames")[14], "CoA(12:0)|PA(14:0/12:0)")
    expect_equal(igraph::E(g_new)$value[14], 0.8)
})

## function add_edge_attributes_from_data.frame
test_that("add_edge_attributes_from_data.frame works", {
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
    
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
     
    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
     
    ## create the adjacency matrix
    reaction_adj <- create_reaction_adjacency_matrix(reaction_l)
     
    ## create graph
    g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)
     
    ## attributes: data.frame
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
            c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
    ))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, attributes, column_from = "from", 
       column_to = "to")

    expect_equal(igraph::E(g_new)$weight, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, 0.8, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$weight[12], 0.5)
    expect_equal(attr(igraph::E(g_new), "vnames")[14], "CoA(12:0)|PA(14:0/12:0)")
    expect_equal(igraph::E(g_new)$weight[14], 0.8)
    expect_error(
        LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, attributes = attributes, 
            column_from = "foo", column_to = "to"), 
        " not in ")
    expect_error(
        LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, attributes = attributes, 
            column_from = "from", column_to = "foo"), 
        " not in ")
    
    ## some vertex names not in g
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
            c("CoA(12:0)", "foo", 0.8),
            c("foo", "PA(14:0/12:0)", 0.8)))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, 
        attributes, column_from = "from", column_to = "to")
    expect_equal(igraph::E(g_new)$weight, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$weight[12], 0.5)
    
    ## all vertex names not in g
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "foo", 0.8),
            c("foo", "PA(14:0/12:0)", 0.8)))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, 
        attributes, column_from = "from", column_to = "to")
    expect_equal(igraph::E(g_new)$weight, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)))
})

## function add_edge_attributes
test_that("add_edge_attributes works", {
    
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
     
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
     
    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
     
    ## create the adjacency matrix
    reaction_adj <- create_reaction_adjacency_matrix(reaction_l)
    
    ## create graph
    g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)

    ## attributes: data.frame
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
            c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
    ))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes, 
        cols_vertex = c("from", "to"))
    expect_equal(igraph::E(g_new)$weight, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, 0.8, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$weight[12], 0.5)
    expect_equal(attr(igraph::E(g_new), "vnames")[14], "CoA(12:0)|PA(14:0/12:0)")
    expect_equal(igraph::E(g_new)$weight[14], 0.8)
    expect_error(
        add_edge_attributes(g = g, attributes = attributes, 
            cols_vertex = "name"), 
        "'cols_vertex' has to be of length 2.")
    
    ## some vertex names not in g
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
            c("CoA(12:0)", "foo", 0.8),
            c("foo", "PA(14:0/12:0)", 0.8)))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes, 
        cols_vertex = c("from", "to"))
    expect_equal(igraph::E(g_new)$weight, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$weight[12], 0.5)
    
    ## all vertex names not in g
    attributes <- data.frame(
        rbind(
            c("CoA(12:0)", "foo", 0.8),
            c("foo", "PA(14:0/12:0)", 0.8)))
    names(attributes) <- c("from", "to", "weight")
    attributes$weight <- as.numeric(attributes$weight)
    
    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes)
    expect_equal(igraph::E(g_new)$weight, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)))
    
    ## attributes: matrix
    attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
        dimnames = list(
            c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)"),
            c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)")))
    
    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes)
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$value[12], 0.5)
    expect_equal(attr(igraph::E(g_new), "vnames")[14], "CoA(12:0)|PA(14:0/12:0)")
    expect_equal(igraph::E(g_new)$value[14], 0.8)
    
    ## some vertex names not in g
    attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
        dimnames = list(
            c("CoA(12:0)", "foo", "PA(14:0/12:0)"),
            c("CoA(12:0)", "PA(12:0/0:0)", "foo")))

    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes)
    expect_equal(igraph::E(g_new)$value, 
        c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
           NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(attr(igraph::E(g_new), "vnames")[12], "CoA(12:0)|PA(12:0/0:0)")
    expect_equal(igraph::E(g_new)$value[12], 0.5)
    
    ## all vertex names not in g
    attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
        dimnames = list(
            c("foo_1", "foo_2", "foo_3"),
            c("foo_1", "foo_2", "foo_3")))
    
    ## apply the function
    g_new <- add_edge_attributes(g = g, attributes = attributes)
    expect_equal(igraph::E(g_new)$value, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)))
})

## function add_vertex_attributes
test_that("add_vertex_attributes works", {
    FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
     
    ## create data.frame with reactions and reaction order
    reactions <- rbind(
        c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
        c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
        c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
        c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
    )
    reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
        reaction_formula = reactions[, 3], directed = reactions[, 4])
    reactions$order <- as.numeric(reactions$order)
    reactions$directed <- as.logical(reactions$directed)
    
    ## create the list with reactions
    reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)

    ## create the adjacency matrix
    reaction_adj <- create_reaction_adjacency_matrix(reaction_l)
    
    ## create graph
    g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)
    
    ## attributes: data.frame
    attributes_df <- data.frame(
        name = c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)", "DG(12:0/12:0/0:0)",
            "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/12:0/0:0)",
            "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)", "DG(16:0/12:0/0:0)",
            "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)", "FA(12:0)", "FA(14:0)",
            "FA(16:0)", "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)",
            "PA(12:0/16:0)", "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", 
            "PA(14:0/16:0)", "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)",
            "PA(16:0/16:0)"),
        logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
            -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
            -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
            -0.11),
        logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
            -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
            -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
    
    ## apply the function
    g_new <- add_vertex_attributes(g, attributes = attributes_df, 
        col_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        c(NA, NA, NA, -5.08, 0.75, 5.43, -0.62, 2.35, 1.39, 2.91, 0.26, -4.14,
            0.19, 6.18, 0.78, -1.81, 4.66, -0.10, NA, NA, 2.84, -0.81, -0.81,
            -0.32, 0.17, 2.25, -1.94, 0.80, 4.21, 0.20, -3.29, -0.11, NA, NA))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        c(NA, NA, NA, -2.73, 6.14, 1.98, 0.09, 1.57, 1.77, 3.08, 4.04, -3.01,
            1.22, -4.25, 0.39, 0.53, 3.30, 7.10, NA, NA, 2.81, -0.99, -0.09,
            -8.25, 4.94, -3.54, -7.74, -1.98, 0.73, 2.36, 2.53, -0.62, NA, NA))
    expect_error(
        add_vertex_attributes(g = g, attributes = attributes_df, 
            col_vertex = "foo"), 
        " not in ")
    
    ## some vertex names not in g
    attributes_df <- data.frame(
        name = c("foo_1", "foo_2", "CoA(16:0)"),
        logFC_cond1 = c(-5.08,  0.75,  5.43),
        logFC_cond2 = c(-2.73,  6.14,  1.98))
    g_new <- add_vertex_attributes(g, attributes = attributes_df, 
        col_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        c(NA, NA, NA, NA, NA, 5.43, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        c(NA, NA, NA, NA, NA, 1.98, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
    
    ## all vertex names not in g
    attributes_df <- data.frame(
        name = c("foo_1", "foo_2", "foo_3"),
        logFC_cond1 = c(-5.08,  0.75,  5.43),
        logFC_cond2 = c(-2.73,  6.14,  1.98))
    g_new <- add_vertex_attributes(g, attributes = attributes_df, 
        col_vertex = "name")
    expect_equal(igraph::V(g_new)$logFC_cond1, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA)))
    expect_equal(igraph::V(g_new)$logFC_cond2, 
        as.numeric(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
            NA, NA, NA, NA)))
    
})
