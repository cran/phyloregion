## -----------------------------------------------------------------------------
library(phyloregion)

## ---- fig.align="left", fig.cap="__Figure 2.__ Phylogenetic tree of the woody plants of southern Africa inferred from DNA barcodes using a maximum likelihood approach and transforming branch lengths to millions of years ago by enforcing a relaxed molecular clock and multiple calibrations [@Daru2015ddi]."----
library(ape)
library(Matrix)
library(terra)
data(africa)
sparse_comm <- africa$comm

tree <- africa$phylo
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
par(mar=c(2,2,2,2))
plot(tree, show.tip.label=FALSE)

## -----------------------------------------------------------------------------
s <- vect(system.file("ex/nigeria.json", package="phyloregion"))

set.seed(1)
m <- as.data.frame(spatSample(s, 1000, method = "random"),
                   geom = "XY")[-1]
names(m) <- c("lon", "lat")
species <- paste0("sp", sample(1:100))
m$taxon <- sample(species, size = nrow(m), replace = TRUE)

pt <- points2comm(dat = m, res = 0.5, lon = "lon", lat = "lat",
            species = "taxon") # This generates a list of two objects
head(pt[[1]][1:5, 1:5])

## -----------------------------------------------------------------------------
s <- vect(system.file("ex/nigeria.json", package="phyloregion"))
sp <- random_species(100, species=5, pol=s)
pol <- polys2comm(dat = sp)
head(pol[[1]][1:5, 1:5])

## -----------------------------------------------------------------------------
fdir <- system.file("NGAplants", package="phyloregion")
files <- file.path(fdir, dir(fdir))
ras <- rast2comm(files) 
head(ras[[1]][1:5, 1:5])

## ---- fig.align="left", fig.cap="__Figure 3.__ Species richness of plants in Nigeria across equal area grid cells. This is to demonstrate how the function `plot` works.", out.width = '50%'----
s <- vect(system.file("ex/SR_Naija.json", package="phyloregion"))
par(mar=rep(0,4))
plot(s, "SR", border=NA, type = "continuous", 
     col = hcl.colors(20, palette = "Blue-Red 3", rev=FALSE))

## ---- eval=TRUE---------------------------------------------------------------
library(Matrix) 
data(africa)
sparse_comm <- africa$comm
dense_comm <- as.matrix(sparse_comm) 
object.size(dense_comm)
object.size(sparse_comm)

## ---- fig.align="left", fig.cap="__Figure 4.__ Geographic distributions of weighted endemism for woody plants of southern Africa.", out.width = '50%'----
library(terra)
data(africa)
p <- vect(system.file("ex/sa.json", package = "phyloregion"))
Endm <- weighted_endemism(africa$comm)
m <- merge(p, data.frame(grids=names(Endm), WE=Endm), by="grids")
m <- m[!is.na(m$WE),]

par(mar=rep(0,4))
plot(m, "WE", col = hcl.colors(20, "Blue-Red 3"), 
     type="continuous", border = NA)

## ---- fig.align="left", fig.cap="__Figure 5.__ Geographic distributions of phylogenetic diversity for woody plants of southern Africa.", out.width = '50%'----
data(africa)
comm <- africa$comm
tree <- africa$phylo
poly <- vect(system.file("ex/sa.json", package = "phyloregion"))

mypd <- PD(comm, tree)
head(mypd)

M <- merge(poly, data.frame(grids=names(mypd), pd=mypd), by="grids")
M <- M[!is.na(M$pd),]
head(M)

par(mar=rep(0,4))
plot(M, "pd", border=NA, type="continuous",
            col = hcl.colors(20, "Blue-Red 3"))

## ---- fig.align="left", fig.cap="__Figure 6.__ Geographic distributions of phylogenetic endemism for woody plants of southern Africa.", out.width = '50%'----
library(terra)
data(africa)
comm <- africa$comm
tree <- africa$phylo
poly <- vect(system.file("ex/sa.json", package = "phyloregion"))

pe <- phylo_endemism(comm, tree)
head(pe)

mx <- merge(poly, data.frame(grids=names(pe), pe=pe), by="grids")
mx <- mx[!is.na(mx$pe),]
head(mx)

par(mar=rep(0,4))
plot(mx, "pe", border=NA, type="continuous",
            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))

## ---- fig.align="left", fig.cap="__Figure 7.__ Geographic distributions of evolutionary distinctiveness and global endangerment for woody plants of southern Africa.", out.width = '50%'----
data(africa)
comm <- africa$comm
threat <- africa$IUCN
tree <- africa$phylo
poly <- vect(system.file("ex/sa.json", package = "phyloregion"))

x <- EDGE(threat, tree, Redlist = "IUCN", species="Species")
head(x)

y <- map_trait(comm, x, FUN = sd, pol=poly)

par(mar=rep(0,4))
plot(y, "traits", border=NA, type="continuous",
            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))

## -----------------------------------------------------------------------------
data(africa)
p <- vect(system.file("ex/sa.json", package = "phyloregion"))
sparse_comm <- africa$comm

tree <- africa$phylo
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
pb <- phylobeta(sparse_comm, tree)

## ---- message=FALSE, results='hide', warning=FALSE----------------------------
y <- phyloregion(pb[[1]], pol=p)

## -----------------------------------------------------------------------------
plot_NMDS(y, cex=3)
text_NMDS(y)

par(mar=rep(0,4))
plot(y, palette="NMDS")

## ---- eval=TRUE---------------------------------------------------------------
sessionInfo()

