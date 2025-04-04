% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_reaction.R
\name{create_reaction_adjacency_matrix}
\alias{create_reaction_adjacency_matrix}
\title{Create adjacency matrix from reactions}
\usage{
create_reaction_adjacency_matrix(reaction_l)
}
\arguments{
\item{reaction_l}{list as obtained from \code{create_reactions}}
}
\value{
matrix
}
\description{
The function \code{create_reaction_adjacency_matrix} creates an adjacency
matrix connecting substrates and products that are linked by reactions.
}
\details{
The function \code{create_reaction_adjacency_matrix} accepts the output of 
the \code{create_reactions} function.

The adjacency matrix can be used in subsequent analysis for network analysis,
e.g. by converting the adjacency matrix to a graph via
\code{igraph::graph_from_adjacency_matrix}.
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
create_reaction_adjacency_matrix(reaction_l = reaction_l)
}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}
}
