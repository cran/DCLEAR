% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.r
\name{simulate_core}
\alias{simulate_core}
\title{simulate_core}
\usage{
simulate_core(config, tree, mutation_site, outcome)
}
\arguments{
\item{config}{simulation configuration; a lineage_tree_config object}

\item{tree}{a matrix representing the lineage tree}

\item{mutation_site}{a binary matrix indicating the mutation sites}

\item{outcome}{a character matrix}
}
\value{
a `lineage_tree` object
}
\description{
Simulate a cell lineage tree
Adoped from https://github.com/elifesciences-publications/CRISPR_recorders_sims/blob/master/MATLAB_sims/GESTALT_30hr_1x_simulation.m
}
