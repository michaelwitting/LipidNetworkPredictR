## function collapse_network
test_that("collapse_network works", {
    
    ## undirected graph
    ## create an example graph
    edges <- data.frame(from = c("A", "B", "C", "D", "C", "H"), 
        to = c("B", "C", "D", "E", "H", "I"))
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    plot(g, vertex.size = 30, edge.width = 5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E", "I")  
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D", "H", "E", "I"))
    expect_equal(names(V(collapsed_g)), c("A", "D", "I", "E"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(7, 7))
    expect_equal(rownames(adj), c("A", "B", "C", "D", "H", "E", "I"))
    expect_equal(as.vector(adj[, 1]), c(0, 1, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 1, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 1, 0, 0, 1, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 1, 0, 0, 0, 1))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 0, 1, 0, 0, 0))
    expect_equal(as.vector(adj[, 7]), c(0, 0, 0, 0, 1, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "D", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 1, 1, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 1, 1))
    expect_equal(as.vector(adj[, 3]), c(1, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 1, 0, 0))
    
    ## directed graph (example 1)
    ## create an example graph
    edges <- data.frame(from = c("A", "B", "C", "D", "H", "I"), 
                          to = c("B", "C", "D", "E", "C", "H"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E", "I")  
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(names(V(collapsed_g)), c("A", "D", "E", "I"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(7, 7))
    expect_equal(rownames(adj), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 0, 0, 0, 1, 0))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 7]), c(0, 0, 0, 1, 0, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "D", "E", "I"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 1))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 0, 0))
    
    ## directed graph (example 2)
    ## create an example graph
    edges <- data.frame(from = c("A", "B", "C", "D", "C", "H"), 
                          to = c("B", "C", "D", "E", "H", "I"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E", "I")  
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D", "H", "E", "I"))
    expect_equal(names(V(collapsed_g)), c("A", "D", "I", "E"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(7, 7))
    expect_equal(rownames(adj), c("A", "B", "C", "D", "H", "E", "I"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 0, 1, 0, 0, 0))
    expect_equal(as.vector(adj[, 7]), c(0, 0, 0, 0, 1, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "D", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 0, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 1, 0, 0))
    
    ## directed graph (example 3)
    ## create an example graph
    edges <- data.frame(from = c("A", "B", "C", "D", "C", "H", "H", "I"), 
                          to = c("B", "C", "D", "E", "H", "C", "I", "H"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E", "I")  
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(names(V(collapsed_g)), c("A", "D", "I", "E"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(7, 7))
    expect_equal(rownames(adj), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 1, 0, 0, 1, 0))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 7]), c(0, 0, 0, 1, 0, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "D", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 1, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 0, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 1, 0, 0))
    
    ## directed graph (example 4)
    ## create an example graph
    edges <- data.frame(from = c("A", "B", "C", "D", "C", "H", "H", "I"), 
                          to = c("B", "C", "D", "E", "H", "C", "I", "H"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E", "H", "I")  
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(names(V(collapsed_g)), c("A", "D", "H", "E", "I"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(7, 7))
    expect_equal(rownames(adj), c("A", "B", "C", "D", "H", "I", "E"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(0, 0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 1, 0, 0, 1, 0))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 7]), c(0, 0, 0, 1, 0, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(5, 5))
    expect_equal(rownames(adj), c("A", "D", "H", "E", "I"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 0, 0, 0, 1))
    expect_equal(as.vector(adj[, 4]), c(0, 1, 0, 0, 0))
    
    ## directed graph (example 5)
    ## create an example graph
    edges <- data.frame(from = c("A", "A", "B", "C"), 
                          to = c("B", "D", "C", "D"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "C")
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D"))
    expect_equal(names(V(collapsed_g)), c("A", "C"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "B", "C", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(1, 0, 1, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(2, 2))
    expect_equal(rownames(adj), c("A", "C"))
    expect_equal(as.vector(adj[, 1]), c(0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0))
    
    ## directed graph (example 6)
    ## create an example graph
    edges <- data.frame(from = c("A", "A", "B", "C"), 
                          to = c("B", "D", "C", "D"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "C", "D")
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "D"))
    expect_equal(names(V(collapsed_g)), c("A", "C", "D"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(4, 4))
    expect_equal(rownames(adj), c("A", "B", "C", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(1, 0, 1, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(3, 3))
    expect_equal(rownames(adj), c("A", "C", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 1, 0))
    
    ## directed graph (example 7)
    ## create an example graph
    edges <- data.frame(from = c("A", "A", "B", "C", "E", "F"), 
                          to = c("B", "E", "C", "D", "F", "D"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "B")
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "E", "F", "D"))
    expect_equal(names(V(collapsed_g)), c("A", "B", "D"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(6, 6))
    expect_equal(rownames(adj), c("A", "B", "C", "E", "F", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(3, 3))
    expect_equal(rownames(adj), c("A", "B", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 1, 0))
    
    ## directed graph (example 7)
    ## create an example graph
    edges <- data.frame(from = c("A", "A", "B", "C", "E", "F"), 
                          to = c("B", "E", "C", "D", "F", "D"))
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    plot(g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## assign data availability to nodes, TRUE for measured nodes
    igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E")
    
    ## collapse the network
    collapsed_g <- collapse_network(g, .vertex_attr_names = "data_available")
    plot(collapsed_g, vertex.size = 30, edge.width = 2, edge.arrow.size = 1.5)
    
    ## check vertex names
    expect_equal(names(V(g)), c("A", "B", "C", "E", "F", "D"))
    expect_equal(names(V(collapsed_g)), c("A", "E", "D"))
    
    ## check adjacency matrix (edges)
    adj <- igraph::as_adjacency_matrix(graph = g, sparse = FALSE)
    expect_equal(dim(adj), c(6, 6))
    expect_equal(rownames(adj), c("A", "B", "C", "E", "F", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(0, 1, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 4]), c(1, 0, 0, 0, 0, 0))
    expect_equal(as.vector(adj[, 5]), c(0, 0, 0, 1, 0, 0))
    expect_equal(as.vector(adj[, 6]), c(0, 0, 1, 0, 1, 0))
    
    adj <- igraph::as_adjacency_matrix(graph = collapsed_g, sparse = FALSE)
    expect_equal(dim(adj), c(3, 3))
    expect_equal(rownames(adj), c("A", "E", "D"))
    expect_equal(as.vector(adj[, 1]), c(0, 0, 0))
    expect_equal(as.vector(adj[, 2]), c(1, 0, 0))
    expect_equal(as.vector(adj[, 3]), c(1, 1, 0))
})
