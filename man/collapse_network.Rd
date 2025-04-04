% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_analysis.R
\name{collapse_network}
\alias{collapse_network}
\title{Collapse network by eliminating vertices with unmeasured values}
\usage{
collapse_network(g, .vertex_attr_names)
}
\arguments{
\item{g}{\code{igraph} object}

\item{.vertex_attr_names}{character of length 1, defining the 
\code{vertex_attr_names} of \code{graph} to use for specifying if there 
is a measurement available for the given vertex}
}
\description{
The function \code{collapse_network} will remove vertices without data by
collapsing their paths while preserving connections between measured 
vertices.
}
\details{
As an example, in the graph \code{A <-> B <-> C <-> D}. The vertices
\code{B} and \code{C} are not measured. The resulting graph after 
applying the function will be \code{A <-> D}.
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

g <- igraph::graph_from_adjacency_matrix(reaction_adj, weighted = TRUE, 
    diag = FALSE)

## attribute_type: vertex
attributes_df <- data.frame(
   name = c("CoA (12:0", "CoA (14:0", "CoA (16:0", "DG 12:0/12:0/0:0",
       "DG 12:0/14:0/0:0", "DG 12:0/16:0/0:0", "DG 14:0/12:0/0:0",
       "DG 14:0/14:0/0:0", "DG 14:0/16:0/0:0", "DG 16:0/12:0/0:0",
       "DG 16:0/14:0/0:0", "DG 16:0/16:0/0:0", "FA (12:0", "F (14:0",
       "FA (16:0", "PA (12:0/0:0", "PA 12:0/12:0", "PA 12:0/14:0",
       "PA 12:0/16:0", "PA 14:0/0:0", "PA 14:0/12:0", "PA 14:0/14:0", 
       "PA 14:0/16:0", "PA 16:0/0:0", "PA 16:0/12:0", "PA 16:0/14:0",
       "PA 16:0/16:0"),
   logFC_cond1 = c(-5.08,  0.75,  5.43, -0.62,  2.35, 1.39, 2.91,  0.26, 
       -4.14,  0.19,  6.18, 0.78, -1.81,  4.66, -0.10,  2.84, -0.81,
       -0.81, -0.32,  0.17,  2.25, -1.94,  0.80, 4.21,  0.20, -3.29, 
       -0.11),
   logFC_cond2 = c(-2.73,  6.14,  1.98,  0.09,  1.57,  1.77,  3.08,  4.04,
       -3.01, 1.22, -4.25, 0.39, 0.53, 3.30, 7.10, 2.81, -0.99, -0.09,
       -8.25, 4.94, -3.54, -7.74, -1.98, 0.73,  2.36,  2.53, -0.62))
       
## add vertex attributes to g
g <- add_attributes(g, attribute_type = "vertex", attributes = attributes_df, 
    cols_vertex = "name")

## collapse the network
g_collapsed <- collapse_network(g = g, .vertex_attr_names = "logFC_cond1")

## plot the original and collapsed graphs
par(mfrow = c(1, 2))
plot(g, main = "Original Graph")
plot(g_collapsed, main = "Collapsed Graph")
}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}
}
