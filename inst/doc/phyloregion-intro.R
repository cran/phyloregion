## -----------------------------------------------------------------------------
library(phyloregion)

## ---- eval=TRUE---------------------------------------------------------------
library(Matrix) 
data(africa)
sparse_comm <- africa$comm
dense_comm <- as.matrix(sparse_comm) 
object.size(dense_comm)
object.size(sparse_comm)

## -----------------------------------------------------------------------------
tree <- africa$phylo
tree <- ape::keep.tip(tree, colnames(sparse_comm))
pb <- phylobeta(sparse_comm, tree)

## ---- eval=TRUE---------------------------------------------------------------
sessionInfo()

