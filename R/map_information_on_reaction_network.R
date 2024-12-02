#' @name add_attributes
#' 
#' @title Add attributes to igraph object
#' 
#' @description 
#' The function \code{add_attributes} adds attributes to an \code{igraph} 
#' object. Attributes can be either added for edges 
#' \code{attribute_type = "edges"} or vertices 
#' \code{attribute_type = "vertices"}. The function will return a \code{igraph}
#' object.
#' 
#' For \code{attribute_type = "edges"}, the function adds edge weights to a 
#' \code{igraph} object \code{g}. The weights are stored in the 
#' \code{attributes} object. The function will return a \code{igraph} object 
#' with updated edge weights. 
#' The \code{attributes} object can be either a \code{matrix} or a
#' \code{data.frame}. 
#' The \code{matrix} is an adjacency matrix containing as
#' entries the weights. The weights will be stored in the \code{E(g)$value}
#' slot of the returned \code{igraph} object.
#' The \code{data.frame} contains the columns \code{vertex}, a 
#' \code{character} vector of length 2, specifying the out- and ingoing 
#' vertices for the edge and the edge weights in the remaining columns. The
#' weights will be stored in the respective slots with same names as the 
#' \code{colnames} of \code{attributes} of the returned \code{igraph} object.
#'
#' For \code{attribute_type = "vertices"}, the function adds vertex attributes 
#' to a \code{igraph} object \code{g}. The values are stored in the 
#' \code{attributes} object. The function will return a \code{igraph} object 
#' with updated vertex attributes. 
#' The \code{attributes} object is a \code{data.frame}. 
#' The \code{data.frame} contains the columns \code{col_vertex}, a 
#' \code{character} vector of length 1, specifying the vertices and
#' the vertex attributes weights in the remaining columns. The
#' attributes will be stored in the respective slots with same names as the 
#' \code{colnames} of \code{attributes} for each vertex 
#' of the returned \code{igraph} object.

#' @details
#' For \code{attribute_type = "edges"}, \code{cols_vertex} has to be adjusted 
#' only when \code{attributes} is a \code{data.frame}. The \code{character} 
#' of length 2 will specify the columns
#' containing the out- and ingoing vertices of the graph.

#' For \code{attribute_type = "vertices"}, \code{cols_vertex} has to be 
#' adjusted such that it specifies the vertices. 
#' The \code{character} of length 1 will specify the column
#' containing the vertices of the graph.
#' 
#' @param g igraph object
#' @param attribute_type \code{character} of length 1, either \code{"edges"} or 
#' \code{"vertices"} 
#' @param attributes \code{data.frame} or \code{matrix} containing edge or 
#' vertex attribute information
#' @param cols_vertex \code{character} of length 1 or length 2 specifying the
#' columns in \code{attributes}
#' 
#' @export
#' 
#' @return igraph object
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @examples 
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
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
#'    name = c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)", "DG(12:0/12:0/0:0)",
#'        "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/12:0/0:0)",
#'        "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)", "DG(16:0/12:0/0:0)",
#'        "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)", "FA(12:0)", "FA(14:0)",
#'        "FA(16:0)", "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)",
#'        "PA(12:0/16:0)", "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", 
#'        "PA(14:0/16:0)", "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)",
#'        "PA(16:0/16:0)"),
#'    logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
#'        -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
#'        -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
#'        -0.11),
#'    logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
#'        -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
#'        -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
#'        
#' ## apply the function
#' add_attributes(g, attribute_type = "vertex", attributes = attributes_df, 
#'     cols_vertex = "name")
#'
#' ## attribute_type: edges, attributes: data.frame
#' attributes <- data.frame(
#'     rbind(
#'         c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
#'         c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
#' ))
#' names(attributes) <- c("from", "to", "weight")
#' 
#' ## apply the function
#' add_attributes(g, attribute_type = "edges", attributes = attributes, cols_vertex = c("from", "to"))
#' 
#' ## attribute_type:edges, attributes: matrix
#' attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
#'     dimnames = list(
#'         c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)"),
#'         c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)")))
#' 
#' ## apply the function
#' add_attributes(g, attribute_type = "edges", attributes = attributes)
#' 
#' 
add_attributes <- function(g, attribute_type = c("edges", "vertex"), 
        attributes, cols_vertex = colnames(attributes)[1:2]) {
    
    ## information can be either submitted via edges or nodes
    ## edges: e.g. correlation between two features
    ## nodes: e.g. log-fold change
    attribute_type <- match.arg(attribute_type)
    
    ## for attribute_type == "edges"
    if (attribute_type == "edges") 
        g <- add_edge_attributes(g = g, attributes = attributes, 
            cols_vertex = cols_vertex)
    
    ## for attribute_type == "vertex"
    if (attribute_type == "vertex")
        g <- add_vertex_attributes(g = g, attributes = attributes, 
            col_vertex = cols_vertex)
        
    ## return the graph
    g
}

#' @name add_edge_attributes_from_data.frame
#' 
#' @title Add edge attributes to \code{igraph} object from \code{data.frame}
#' 
#' @description 
#' The function adds edge weights to a \code{igraph} object \code{g}. The 
#' weights are stored in the \code{attributes} object. The function will return
#' a \code{igraph} object with updated edge weights. 
#' 
#' The \code{attributes} object is a \code{data.frame}. 
#' 
#' The \code{data.frame} contains the columns \code{vertex}, a 
#' \code{character} vector of length 2, specifying the out- and ingoing 
#' vertices for the edge and the edge weights in the remaining columns. The
#' weights will be stored in the respective slots with same names as the 
#' \code{colnames} of \code{attributes} of the returned \code{igraph} object.#' 
#' 
#' @details
#' \code{vertex} has to be adjusted to the \code{attributes} object. 
#' The \code{character} of length 2 will specify the columns
#' containing the out- and ingoing vertices of the graph.
#' 
#' @param g \code{igraph} object
#' @param attributes \code{data.frame} containing edge
#' attribute information
#' @param column_from \code{chararacter} of length 1, specifying the columns 
#' containing the outgoing vertices in \code{attributes} of type
#' \code{data.frame}
#' @param column_to \code{chararacter} of length 1, specifying the columns 
#' containing the ingoing vertices in \code{attributes} of type
#' \code{data.frame}
#' 
#' @return igraph object
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @importFrom stats setNames
#' @importFrom igraph set_edge_attr
#' 
#' @examples 
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
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
#' ## create graph
#' g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)
#' 
#' ## attributes: data.frame
#' attributes <- data.frame(
#'     rbind(
#'         c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
#'         c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
#' ))
#' names(attributes) <- c("from", "to", "weight")
#' 
#' ## apply the function
#' LipidNetworkPredictR:::add_edge_attributes_from_data.frame(g, attributes, column_from = "from", 
#'     column_to = "to")
add_edge_attributes_from_data.frame <- function(g, attributes, 
    column_from = "from", column_to = "to") {
    
    if (!is.data.frame(attributes))
        stop("'attributes' has to be a data.frame.")
        
    if (ncol(attributes) < 3)
        stop("The data.frame must have at least 3 columns: 'from', 'to' and attribute column(s).")
        
    ## extract vertex names and attribute names
    from_vertices <- attributes[[column_from]]
    to_vertices <- attributes[[column_to]]
    cols <- colnames(attributes)
    attribute_names <- cols[!(cols %in% c(column_from, column_to))]
        
    ## iterate through each attribute column to update weights
    for (attr in attribute_names) {
        
        ## create a names vector of new attribute values
        new_values <- stats::setNames(attributes[[attr]], 
            paste(from_vertices, to_vertices, sep = "|"))
            
        ## update the graph's edge attributes
        ## initialize the attribute of the graph
        g <- igraph::set_edge_attr(g, attr, value = NA)
            
        ## match and update values
        g <- igraph::set_edge_attr(g, attr, index = names(new_values), 
            value = as.vector(new_values))
    }
    
    ## return the updated graph
    g
}

#' @name add_edge_attributes
#' 
#' @title Add edge attributes to \code{igraph} object
#' 
#' @description 
#' The function adds edge weights to a \code{igraph} object \code{g}. The 
#' weights are stored in the \code{attributes} object. The function will return
#' a \code{igraph} object with updated edge weights. 
#' 
#' The \code{attributes} object can be either a \code{matrix} or a
#' \code{data.frame}. 
#' 
#' The \code{matrix} is an adjacency matrix containing as
#' entries the weights. The weights will be stored in the \code{E(g)$value}
#' slot of the returned \code{igraph} object.
#' 
#' The \code{data.frame} contains the columns \code{vertex}, a 
#' \code{character} vector of length 2, specifying the out- and ingoing 
#' vertices for the edge and the edge weights in the remaining columns. The
#' weights will be stored in the respective slots with same names as the 
#' \code{colnames} of \code{attributes} of the returned \code{igraph} object.
#' 
#' @details
#' \code{cols_vertex} has to be adjusted only when \code{attributes} is a
#' \code{data.frame}. The \code{character} of length 2 will specify the columns
#' containing the out- and ingoing vertices of the graph.
#' 
#' @param g \code{igraph} object
#' @param attributes \code{matrix} or \code{data.frame} containing edge
#' attribute information
#' @param cols_vertex \code{chararacter} of length 2, specifying the columns 
#' containing the out- and ingoing vertices in \code{attributes} of type
#' \code{data.frame}
#' 
#' @export
#' 
#' @return igraph object
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @importFrom tidyr pivot_longer
#' 
#' @examples
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
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
#' ## create graph
#' g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)
#' 
#' ## attributes: data.frame
#' attributes <- data.frame(
#'     rbind(
#'         c("CoA(12:0)", "PA(12:0/0:0)", 0.5),
#'         c("CoA(12:0)", "PA(14:0/12:0)", 0.8)
#' ))
#' names(attributes) <- c("from", "to", "weight")
#' 
#' ## apply the function
#' add_edge_attributes(g, attributes, cols_vertex = c("from", "to"))
#' 
#' ## attributes: matrix
#' attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
#'     dimnames = list(
#'         c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)"),
#'         c("CoA(12:0)", "PA(12:0/0:0)", "PA(14:0/12:0)")))
#' 
#' ## apply the function
#' add_edge_attributes(g, attributes)
add_edge_attributes <- function(g, attributes, cols_vertex = colnames(attributes)[1:2]) {
    
    ## information can be submitted via edges
    ## edges: e.g. correlation between two features
    
    ## there are two ways of storing the edge information: either in a long
    ## data.frame or in an adjacency matrix
    ## long data.frame: the columns values_from and values_to contain the edge 
    ## information while the other columns (at least one) contain the weights
    if (is.data.frame(attributes)) {
        if (length(cols_vertex) != 2)
            stop("'cols_vertex' has to be of length 2.")
        
        if (!all(cols_vertex %in% colnames(attributes))) 
            stop("'cols_vertex' have to be in colnames(attributes).")
            
        g <- add_edge_attributes_from_data.frame(g, attributes = attributes,
            column_from = cols_vertex[1], column_to = cols_vertex[2])
    }
    
    ## adjacency matrix: the matrix will be first converted into a long 
    ## data.frame, the graph is updated using the long data.frame
    if (is.matrix(attributes)) {
        if (ncol(attributes) != nrow(attributes)) 
            stop("'attributes' has to be a square matrix.")
        
        ## convert into long data.frame
        attributes_df <- attributes |>
            as.data.frame()
        attributes_df[["Row"]] <- rownames(attributes)
        attributes_df <- tidyr::pivot_longer(attributes_df, 
            cols = -c("Row"), names_to = "Col", values_to = "value")
        
        ## add edge attributes from long data.frame
        g <- add_edge_attributes_from_data.frame(g, 
            attributes = attributes_df, 
            column_from = "Row", column_to = "Col")
    }
    
    ## return the updated graph object
    g
}

#' @name add_vertex_attributes
#' 
#' @title Add vertex attributes to \code{igraph} object
#' 
#' @description 
#' The function adds vertex attributes to a \code{igraph} object \code{g}. The 
#' values are stored in the \code{attributes} object. The function will return
#' a \code{igraph} object with updated vertex attributes. 
#' 
#' The \code{attributes} object is a \code{data.frame}. 
#' 
#' The \code{data.frame} contains the columns \code{col_vertex}, a 
#' \code{character} vector of length 1, specifying the vertices and
#' the vertex attributes weights in the remaining columns. The
#' attributes will be stored in the respective slots with same names as the 
#' \code{colnames} of \code{attributes} for each vertex 
#' of the returned \code{igraph} object.
#' 
#' @details
#' \code{col_vertex} has to be adjusted such that it specifies the vertices. 
#' The \code{character} of length 1 will specify the column
#' containing the vertices of the graph.
#' 
#' @param g \code{igraph} object
#' @param attributes \code{data.frame} containing vertex
#' attribute information
#' @param col_vertex \code{chararacter} of length 1, specifying the column 
#' containing  vertices in \code{attributes} of type \code{data.frame}
#' 
#' @export
#' 
#' @return igraph object
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @importFrom igraph set_vertex_attr
#' 
#' @examples
#' FA <- c("FA(12:0)", "FA(14:0)", "FA(16:0)")
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
#' ## create graph
#' g <- igraph::graph_from_adjacency_matrix(reaction_adj, mode = "directed", weighted = TRUE)
#' 
#' ## attributes: data.frame
#' attributes_df <- data.frame(
#'    name = c("CoA(12:0)", "CoA(14:0)", "CoA(16:0)", "DG(12:0/12:0/0:0)",
#'        "DG(12:0/14:0/0:0)", "DG(12:0/16:0/0:0)", "DG(14:0/12:0/0:0)",
#'        "DG(14:0/14:0/0:0)", "DG(14:0/16:0/0:0)", "DG(16:0/12:0/0:0)",
#'        "DG(16:0/14:0/0:0)", "DG(16:0/16:0/0:0)", "FA(12:0)", "FA(14:0)",
#'        "FA(16:0)", "PA(12:0/0:0)", "PA(12:0/12:0)", "PA(12:0/14:0)",
#'        "PA(12:0/16:0)", "PA(14:0/0:0)", "PA(14:0/12:0)", "PA(14:0/14:0)", 
#'        "PA(14:0/16:0)", "PA(16:0/0:0)", "PA(16:0/12:0)", "PA(16:0/14:0)",
#'        "PA(16:0/16:0)"),
#'    logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
#'        -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
#'        -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
#'        -0.11),
#'    logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
#'        -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
#'        -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
#' 
#' ## apply the function
#' add_vertex_attributes(g = g, attributes = attributes_df, col_vertex = "name")
add_vertex_attributes <- function(g, attributes, col_vertex = colnames(attributes_df)[1]) {
    
    ## check arguments
    if (!is.data.frame(attributes))
        stop("'attributes' has to be a data.frame.")
    
    if (length(col_vertex) != 1)
        stop("'col_vertex' has to be of length 1.")
    
    if (any(duplicated(attributes[[col_vertex]])))
        stop("'attributes[[col_vertex]]' contains duplicated entries.")
    
    ## obtain the colnames, obtain the columns that store attribute information
    cols_df <- colnames(attributes)
    cols_df_attr <- cols_df[cols_df != col_vertex]
    
    ## iterate trough the attribute columns and add attribute information
    for (attr_i in cols_df_attr) {
        ## set V(g) for cols_df_attr to NA
        g <- igraph::set_vertex_attr(graph = g, name = attr_i, value = NA)
        ## set V(g) for cols_df_attr to attribute values
        g <- igraph::set_vertex_attr(graph = g, name = attr_i, 
            index = attributes[[col_vertex]], 
            value = attributes[[attr_i]])   
    }
    
    ## return the graph
    g
    
}
