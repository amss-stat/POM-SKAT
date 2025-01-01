
# POM-SKAT Implementation

## Overview

This script implements the POM-SKAT (Proportional Odds Model based Sequence Kernel Association Test) to test the association between genetic variants and phenotypes.

## Input

Genotype Matrix (G)

Dimensions: n $\times$ m

Represents n individuals and m genetic variants.

Phenotype Vector (Y)

Dimensions: n

Contains the phenotypic outcomes for n individuals.

Covariate Matrix (X)

Dimensions: n $\times$ p

Includes p covariates for n individuals.

## Dependencies

Ensure the following R packages are installed:

MASS;PearsonDS

## Outputs

P-value (p): Significance of the association.
