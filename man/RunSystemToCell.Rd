% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunSystemToCell.R
\name{RunSystemToCell}
\alias{RunSystemToCell}
\title{RunSystemToCell}
\usage{
RunSystemToCell(
  sys.small,
  ground.truth,
  assay,
  blend = "sum",
  meta.data.to.map
)
}
\arguments{
\item{sys.small}{A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.}

\item{ground.truth}{Ground truth signaling mechanisms present in sys.small.}

\item{assay}{The assay to run the SystemToCell transformation on. Defaults to "RNA."}

\item{blend}{Choice of linear operator to combine edges. Defaults to "sum", also accepts "mean"}

\item{meta.data.to.map}{A character vector of metadata names present in the original object which will be carried to the NICHES objects}
}
\description{
Condenses signaling edges landing on each cell within a Seurat object, from all other cells in the system treated as a single source. 
Outputs another Seurat object, but where the rows of the matrix are ligand-receptor mechanisms
and the columns are each a single receiving cell barcode. The information in the matrix is a sum (or an average, depending on user preference) of
all signaling edges landing on that particular cell, from all cells in the system (including from itself.)
This transformation allows rapid manipulation and dimensional reduction of how a cell is connected within the system.
The default assay of this object is called "SystemToCell" to distinguish it from other Seurat objects.
Meta.data slots by default contain "ReceivingType" information, which is the celltypes for each point, 
and "ReceivingCell" which is the exact cell barcode present in the original Seurat object.
}
