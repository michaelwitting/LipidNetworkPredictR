#' @name get_valid_edges
#' 
#' @title Get valid edges
#' 
#' @description
#' The function \code{get_valid_edges} gets valid edges from a graph.
#' A valid edge is defined as follows:
#' (1) take from a graph g the two measured vertices with indices i, j; 
#' (2) obtain all shortest paths between the measured vertices;
#' (3) assign edge = 1 (valid), if there is/are ONLY one/multiple unmeasured  
#' intermediate vertex/vertices or ZERO intermediate vertices within the shortest
#' paths
#' (4) do not assign edge = 1 (not valid), if there is AT LEAST one measured
#' intermediate vertex within shortest path
#' 
#' @details 
#' Measured vertices will be denoted via \code{vertex_attr_logical}.
#' 
#' The function \code{get_valid_edges} is a helper function for 
#' \code{collapse_network}.
#' 
#'  
#' @param vertices vertices of \code{igraph} object
#' @param paths list of \code{igraph.vs} objects
#' @param vertex_attr_logical logical vector storing information if a vertex is
#' measured
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @examples
#' edges <- data.frame(from = c("A", "A", "B", "C", "E", "F"), 
#'     to = c("B", "E", "C", "D", "F", "D"))
#' g <- igraph::graph_from_data_frame(edges, directed = TRUE)
#' ## assign data availability to nodes, TRUE for measured vertices
#' igraph::V(g)$data_available <- igraph::V(g)$name %in% c("A", "D", "E")
#' 
#' ## check if graph is directed
#' is_directed <- igraph::is_directed(graph = g)
#' 
#' ## obtain the vertex attributes and create a vector, vertex_attr_logical,
#' ## that stores information if a measurement is available (TRUE) or
#' ## not (FALSE)
#' vertex_attr <- igraph::vertex_attr(graph = g, name = "data_available")
#' vertex_attr_logical <- logical(length(vertex_attr))
#' 
#' if (is.numeric(vertex_attr) | is.character(vertex_attr))
#'     vertex_attr_logical[!is.na(vertex_attr)] <- TRUE
#' if (is.logical(vertex_attr))
#'     vertex_attr_logical[vertex_attr] <- TRUE
#' 
#' ## obtain the vertices where information is available
#' vertices_measured <- igraph::V(g)[vertex_attr_logical]
#' 
#' vertex_1 <- vertices_measured[1]
#' vertex_2 <- vertices_measured[3]
#' 
#' ## find the shortest path from vertex_1 to vertex_2
#' path_1 <- suppressWarnings(igraph::all_shortest_paths(graph = g, 
#'     from = vertex_1, to = vertex_2)$res)
#' 
#' ## find the shortest path from vertex_2 to vertex_1
#' path_2 <- suppressWarnings(igraph::all_shortest_paths(graph = g, 
#'     from = vertex_2, to = vertex_1)$res)
#' 
#' ## check and add paths if valid, add a valid path when between 
#' ## shortest paths there are only "unmeasured" vertices
#' ## apply the helper function add_valid_edges for path_1 and path_2
#' get_valid_edges(vertices = igraph::V(g), paths = path_1, 
#'     vertex_attr_logical = vertex_attr_logical)
#' get_valid_edges(vertices = igraph::V(g), paths = path_2,
#'     vertex_attr_logical = vertex_attr_logical)
get_valid_edges <- function(vertices, paths, vertex_attr_logical) {
    lapply(paths, function(path) {
        len_path <- length(path)
        if (len_path > 1 && 
            all(!vertex_attr_logical[match(path, vertices)][-c(1, len_path)])) {
            c(vertices[path[1]]$name, vertices[path[len_path]]$name)
        }
    })
}

#' @name collapse_network
#' 
#' @title Collapse network by eliminating vertices with unmeasured values
#' 
#' @description
#' The function \code{collapse_network} will remove vertices without data by
#' collapsing their paths while preserving connections between measured 
#' vertices.
#' 
#' @details
#' As an example, in the graph \code{A <-> B <-> C <-> D}. The vertices
#' \code{B} and \code{C} are not measured. The resulting graph after 
#' applying the function will be \code{A <-> D}.
#'  
#' @param g \code{igraph} object
#' @param .vertex_attr_names character of length 1, defining the 
#' \code{vertex_attr_names} of \code{graph} to use for specifying if there 
#' is a measurement available for the given vertex
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom methods is
#' @importFrom utils combn
#' @importFrom igraph vertex_attr_names V all_shortest_paths graph_from_edgelist
#' @importFrom igraph is_directed
#' 
#' @examples
#' FA <- c("FA 12:0", "FA 14:0", "FA 16:0")
#' 
#' ## create data.frame with reactions and reaction order
#' reactions <- rbind(
#'     c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
#'     c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
#'     c(3, "RHEA:19709", "M_AcylCoA + M_LPA = M_CoA + M_PA", FALSE),
#'     c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
#' )
#' reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
#'     reaction_formula = reactions[, 3], directed = reactions[, 4])
#' reactions$order <- as.numeric(reactions$order)
#' reactions$directed <- as.logical(reactions$directed)
#' 
#' ## create the list with reactions
#' reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)
#' 
#' ## create the adjacency matrix
#' reaction_adj <- create_reaction_adjacency_matrix(reaction_l)
#' 
#' g <- igraph::graph_from_adjacency_matrix(reaction_adj, weighted = TRUE, 
#'     diag = FALSE)
#' 
#' ## attribute_type: vertex
#' attributes_df <- data.frame(
#'    name = c("CoA (12:0", "CoA (14:0", "CoA (16:0", "DG 12:0/12:0/0:0",
#'        "DG 12:0/14:0/0:0", "DG 12:0/16:0/0:0", "DG 14:0/12:0/0:0",
#'        "DG 14:0/14:0/0:0", "DG 14:0/16:0/0:0", "DG 16:0/12:0/0:0",
#'        "DG 16:0/14:0/0:0", "DG 16:0/16:0/0:0", "FA (12:0", "F (14:0",
#'        "FA (16:0", "PA (12:0/0:0", "PA 12:0/12:0", "PA 12:0/14:0",
#'        "PA 12:0/16:0", "PA 14:0/0:0", "PA 14:0/12:0", "PA 14:0/14:0", 
#'        "PA 14:0/16:0", "PA 16:0/0:0", "PA 16:0/12:0", "PA 16:0/14:0",
#'        "PA 16:0/16:0"),
#'    logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
#'        -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
#'        -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
#'        -0.11),
#'    logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
#'        -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
#'        -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
#'        
#' ## add vertex attributes to g
#' g <- add_attributes(g, attribute_type = "vertex", attributes = attributes_df, 
#'     cols_vertex = "name")
#'
#' ## collapse the network
#' g_collapsed <- collapse_network(g = g, .vertex_attr_names = "logFC_cond1")
#'
#' ## plot the original and collapsed graphs
#' par(mfrow = c(1, 2))
#' plot(g, main = "Original Graph")
#' plot(g_collapsed, main = "Collapsed Graph")
collapse_network <- function(g, .vertex_attr_names) {
    
    if (methods::is(g) != "igraph")
        stop("'g' has to be an 'igraph' object.")
    
    if (length(.vertex_attr_names) != 1) 
        stop("'.vertex_attr_names' has to be of length 1.")
    
    if (!(.vertex_attr_names %in% igraph::vertex_attr_names(g)))
        stop("'.vertex_attr_names' has to be in 'vertex_attr_names(g)'.")
    
    ## check if graph is directed
    is_directed <- igraph::is_directed(graph = g)
    
    ## obtain the vertex attributes and create a vector, vertex_attr_logical,
    ## that stores information if a measurement is available (TRUE) or
    ## not (FALSE)
    vertex_attr <- igraph::vertex_attr(graph = g, name = .vertex_attr_names)
    vertex_attr_logical <- logical(length(vertex_attr))
    
    if (is.numeric(vertex_attr) | is.character(vertex_attr))
        vertex_attr_logical[!is.na(vertex_attr)] <- TRUE
    if (is.logical(vertex_attr))
        vertex_attr_logical[vertex_attr] <- TRUE
    
    ## obtain the vertices where information is available
    vertices <- igraph::V(g)
    vertices_measured <- vertices[vertex_attr_logical]
    
    ## generate all pairs of measured vertices
    measured_pairs <- utils::combn(vertices_measured$name, 2, simplify = FALSE)
    
    ## if the network is directed also check in the opposite direction, i.e.
    ## check the path A to B and B to A (if existing)
    if (is_directed) {
        measured_pairs <- c(measured_pairs, 
            lapply(measured_pairs, function(pair) c(pair[2], pair[1])))
    }
    
    ## loop through pairs of measured vertices, find shortest paths and 
    ## validate edges
    new_edges <- lapply(measured_pairs, function(pair) {
        
        vertex_1 <- pair[1]
        vertex_2 <- pair[2]
        
        ## find the shortest path from vertex_1 to vertex_2 / vertex_2 to vertex_1
        paths <- suppressWarnings(
            igraph::all_shortest_paths(g, from = vertex_1, to = vertex_2)$res)
        
        ## check and add paths if valid, add a valid path when between 
        ## shortest paths there are only "unmeasured" vertices
        ## apply the helper function get_valid_edges for path_1 and path_2
        get_valid_edges(paths = paths, vertices = vertices, 
            vertex_attr_logical = vertex_attr_logical)
    })
    
    ## remove NULL values
    new_edges <- unlist(new_edges, recursive = FALSE) |>
        Filter(f = Negate(is.null))
    
    ## if !is_directed the matrix contains e.g. both A -- B and B -- A, remove 
    ## the redundant information
    if (!is_directed) {
        new_edges <- lapply(new_edges, sort)
    }
    
    ## avoid duplicate edges and rbind
    edges_unique <- unique(new_edges) |>
        do.call(what = "rbind")
    
    ## create a simplified/collapsed graph with only measured vertices, 
    ## create a graph from the edgelist
    g_collapsed <- igraph::graph_from_edgelist(el = edges_unique, 
        directed = is_directed)
    
    ## return the object
    g_collapsed
}
