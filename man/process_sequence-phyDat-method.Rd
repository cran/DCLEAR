% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_sequence.r
\name{process_sequence,phyDat-method}
\alias{process_sequence,phyDat-method}
\title{Process sequences}
\usage{
\S4method{process_sequence}{phyDat}(
  x,
  division = 16L,
  dropout_character = "*",
  default_character = "0",
  deletion_character = "-"
)
}
\arguments{
\item{x}{input data in phyDat format}

\item{division}{cell division}

\item{dropout_character}{Dropout character (default: '*')}

\item{default_character}{Default character (default: '0')}

\item{deletion_character}{Deletion character (default: '-')}
}
\value{
a `lineage_tree_config` object
}
\description{
Process sequences
}
\author{
Wuming Gong (gongx030@umn.edu)
}
