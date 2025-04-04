% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_information_on_reaction_network.R
\name{add_edge_attributes}
\alias{add_edge_attributes}
\title{Add edge attributes to \code{igraph} object}
\usage{
add_edge_attributes(g, attributes, cols_vertex = colnames(attributes)[1:2])
}
\arguments{
\item{g}{\code{igraph} object}

\item{attributes}{\code{matrix} or \code{data.frame} containing edge
attribute information}

\item{cols_vertex}{\code{chararacter} of length 2, specifying the columns 
containing the out- and ingoing vertices in \code{attributes} of type
\code{data.frame}}
}
\value{
igraph object
}
\description{
The function adds edge weights to a \code{igraph} object \code{g}. The 
weights are stored in the \code{attributes} object. The function will return
a \code{igraph} object with updated edge weights. 

The \code{attributes} object can be either a \code{matrix} or a
\code{data.frame}. 

The \code{matrix} is an adjacency matrix containing as
entries the weights. The weights will be stored in the \code{E(g)$value}
slot of the returned \code{igraph} object.

The \code{data.frame} contains the columns \code{vertex}, a 
\code{character} vector of length 2, specifying the out- and ingoing 
vertices for the edge and the edge weights in the remaining columns. The
weights will be stored in the respective slots with same names as the 
\code{colnames} of \code{attributes} of the returned \code{igraph} object.
}
\details{
\code{cols_vertex} has to be adjusted only when \code{attributes} is a
\code{data.frame}. The \code{character} of length 2 will specify the columns
containing the out- and ingoing vertices of the graph.

In case of \code{attributes} is a \code{matrix}, the matrix entries will be
stored in \code{E(g)$value}
}
\examples{
FA <- c("FA 12:0", "FA 14:0", "FA 16:0")

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
        c("CoA 12:0", "PA 12:0/0:0", 0.5),
        c("CoA 12:0", "PA 14:0/12:0", 0.8)
))
names(attributes) <- c("from", "to", "weight")
attributes$weight <- as.numeric(attributes$weight)

## apply the function
add_edge_attributes(g, attributes, cols_vertex = c("from", "to"))

## attributes: matrix
attributes <- matrix(c(0, 0.5, 0.8, 0, 0, 0, 0, 0, 0), ncol = 3, byrow = TRUE, 
    dimnames = list(
        c("CoA 12:0", "PA 12:0/0:0", "PA 14:0/12:0"),
        c("CoA 12:0", "PA 12:0/0:0", "PA 14:0/12:0")))

## apply the function
add_edge_attributes(g, attributes)
}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}
}
