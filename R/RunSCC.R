# Direction 1: each runscc can get called, but there will be much redundancy when running RunSCC (e.g. load lr multi times)

# Direction 2: RunSCC wraps the shared functions in each runscc, each runscc only keeps its essential section. User can only interact with RunSCC 
# By doing this, the user can still run specific organization through the parameter switches

# Bug: min.cells.per.ident needs to be passed to each function


#' RunSCC
#' 
#' Performs Single-Cell Connectivity (SCC) transformations on a Seurat object. 
#' By default, references RunCellToCell, RunCellToSystem, and RunSystemToCell functions. 
#' If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
#' Output is a set of specialized Seurat objects (SCC objects) containing cell signaling information.
#' 
#' @param object A Seurat 4.0 object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed. Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 10. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the SCC objects
#' @param position.x The name of the meta.data column specifying location on the spatial x-axis. Only relevant for spatial omics data.
#' @param position.y The name of the meta.data column specifying location on the spatial y-axis. Only relevant for spatial omics data.
#' @param CellToCell 
#' @param CellToSystem 
#' @param SystemToCell 
#' @param CellToCellSpatial 
#' @param CellToNeighborhood 
#' @param NeighborhoodToCell 
#' @param ... Additional parameters to pass to RunCellToCell, RunSystemToCell, RunCellToSystem, or spatial equivalents
#'
#' @export


RunSCC <- function(object,
                        LR.database = 'fantom5',
                        species,
                        assay = 'RNA',
                        min.cells.per.ident = 10,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = T,
                        CellToCellSpatial = T,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = T,
                        rad= 1,
                        ...){
   # Initialize output structure
  output <- list()
  
  # Calculate SCC organizations without spatial restrictions
  
  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(object,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident=min.cells.per.ident,...)}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(object,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident=min.cells.per.ident,...)}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(object,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident=min.cells.per.ident,...)}
  
  # If requested, additionally calculate spatially-limited SCC organizations
  
  if (CellToCellSpatial == T | CellToNeighborhood == T | NeighborhoodToCell == T){
    
    if (is.null(position.x) | is.null(position.y)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    
  ## Define distance between cells here (?), to use for radius calculations. Pass to the below three functions. This will generalize the function to any spatial dataset (?)
  
  }
  
  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(object,assay = assay, species = species,position.x = position.x,position.y = position.y,rad=rad,min.cells.per.ident=min.cells.per.ident,...)} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(object,assay = assay, species = species,position.x = position.x,position.y = position.y,min.cells.per.ident=min.cells.per.ident,...)} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(object,assay = assay, species = species,position.x = position.x,position.y = position.y,min.cells.per.ident=min.cells.per.ident,...)} #Spatially-limited Neighborhood-Cell vectors (niches)

  # Compile objects for output
  return(output)
}

