# LipidNetworkPredictR

## Overview 

Genome scale metabolic networks typically do not represent lipid metabolism in 
full detail. This package aims to use template reactions to generate a lipid 
reaction network at resolution of acyl chains.

## Installation

To install _LipidNetworkPredictR_ from GitHub, install the package via `devtools`:
```r 
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
library(devtools)
install_github("michaelwitting/LipidNetworkPredictR")
```

The code in the Github repository is a place where we publicly disclose our 
development process.

## Quick start

```{r}
library("LipidNetworkPredictR")

## define the fatty acids that the reactions will be built upon
FA <- c("FA(14:0(12Me))", "FA(16:0(14Me))", "FA(15:1(9Z)(14Me))",        
	"FA(17:0(16Me))", "FA(12:0(11Me))", "FA(13:0(12Me))", "FA(14:0(13Me))",
	"FA(15:0(14Me))", "FA(16:0(15Me))", "FA(12:0)", "FA(14:0)")

## create data.frame with reactions and reaction order, for illustrative
## reasons this will only represent a small subset of the actual lipid 
## metabolism
reactions <- rbind(
	c(1, "RHEA:15421", "M_ATP + M_CoA + M_FA = M_PPi + M_AMP + M_AcylCoA", FALSE),
	c(2, "RHEA:15325", "M_Glycerol-3-P + M_AcylCoA = M_CoA + M_LPA", FALSE),
	c(3, "RHEA:19709", "M_LPA + M_AcylCoA = M_CoA + M_PA", FALSE),
	c(4, "RHEA:27429", "M_H2O + M_PA = M_Pi + M_1,2-DG", FALSE)
)
reactions <- data.frame(order = reactions[, 1], reaction_RHEA = reactions[, 2],
	reaction_formula = reactions[, 3], directed = reactions[, 4])
reactions$order <- as.numeric(reactions$order)
reactions$directed <- as.logical(reactions$directed)

## run the function
reaction_l <- create_reactions(substrates = list(FA = FA), reactions = reactions)

## create the adjacency matrix
adj <- create_reaction_adjacency_matrix(reaction_l)
```

The `adj` object can be further analyzed in subsequent analysis (e.g. by using `igraph`).
