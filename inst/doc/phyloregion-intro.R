## -----------------------------------------------------------------------------
library(phyloregion)

## ---- fig.align="left", fig.cap="__Figure 2.__ Phylogenetic tree of the woody plants of southern Africa inferred from DNA barcodes using a maximum likelihood approach and transforming branch lengths to millions of years ago by enforcing a relaxed molecular clock and multiple calibrations [@Daru2015ddi]."----
library(ape)
library(Matrix)
library(sp)
data(africa)
sparse_comm <- africa$comm

tree <- africa$phylo
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
par(mar=c(2,2,2,2))
plot(tree, show.tip.label=FALSE)

## -----------------------------------------------------------------------------
s <- readRDS(system.file("nigeria/nigeria.rds", package = "phyloregion"))

set.seed(1)
m <- data.frame(sp::spsample(s, 10000, type = "nonaligned"))
names(m) <- c("lon", "lat")
species <- paste0("sp", sample(1:1000))
m$taxon <- sample(species, size = nrow(m), replace = TRUE)

pt <- points2comm(dat = m, mask = s, res = 0.5, lon = "lon", lat = "lat",
            species = "taxon")
head(pt[[1]][1:5, 1:5])

## -----------------------------------------------------------------------------
s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
sp <- random_species(100, species=5, shp=s)
pol <- polys2comm(dat = sp, species = "species", trace=0)
head(pol[[1]][1:5, 1:5])

## -----------------------------------------------------------------------------
fdir <- system.file("NGAplants", package="phyloregion")
files <- file.path(fdir, dir(fdir))
ras <- raster2comm(files)
head(ras[[1]])

## ---- fig.align="left", fig.cap="__Figure 3.__ Species richness of plants in Nigeria across equal area grid cells. This is to demonstrate how the function `plot_swatch` works.", out.width = '50%'----
s <- readRDS(system.file("nigeria/SR_Naija.rds", package = "phyloregion"))
par(mar=rep(0,4))
plot_swatch(s, values = s$SR, k = 20, leg=1, border=NA)

## ---- eval=TRUE---------------------------------------------------------------
library(Matrix) 
data(africa)
sparse_comm <- africa$comm
dense_comm <- as.matrix(sparse_comm) 
object.size(dense_comm)
object.size(sparse_comm)

## ---- fig.align="left", fig.cap="__Figure 4.__ Geographic distributions of weighted endemism for woody plants of southern Africa.", out.width = '50%'----
library(raster)
data(africa)
Endm <- weighted_endemism(africa$comm)
head(Endm)
m <- merge(africa$polys, data.frame(grids=names(Endm), WE=Endm), by="grids")
m <- m[!is.na(m@data$WE),]

par(mar=rep(0,4))
plot_swatch(m, values = m$WE, k=20, leg = 3, border = NA)

## ---- fig.align="left", fig.cap="__Figure 5.__ Geographic distributions of phylogenetic diversity for woody plants of southern Africa.", out.width = '50%'----
data(africa)
comm <- africa$comm
tree <- africa$phylo
poly <- africa$polys

mypd <- PD(comm, tree)
head(mypd)

M <- merge(poly, data.frame(grids=names(mypd), pd=mypd), by="grids")
M <- M[!is.na(M@data$pd),]
head(M)

par(mar=rep(0,4))
plot_swatch(M, values = M$pd, k=20, border=NA, leg=3)

## ---- fig.align="left", fig.cap="__Figure 6.__ Geographic distributions of phylogenetic endemism for woody plants of southern Africa.", out.width = '50%'----
library(raster)
data(africa)
comm <- africa$comm
tree <- africa$phylo
poly <- africa$polys

pe <- phylo_endemism(comm, tree)
head(pe)

mx <- merge(poly, data.frame(grids=names(pe), pe=pe), by="grids")
mx <- mx[!is.na(mx@data$pe),]
head(mx)

par(mar=rep(0,4))
plot_swatch(mx, values = mx$pe, k=20, border=NA, leg=3)

## ---- fig.align="left", fig.cap="__Figure 7.__ Geographic distributions of evolutionary distinctiveness and global endangerment for woody plants of southern Africa.", out.width = '50%'----
data(africa)
comm <- africa$comm
threat <- africa$IUCN
tree <- africa$phylo
poly <- africa$polys

x <- EDGE(threat, tree, Redlist = "IUCN", species="Species")
head(x)

y <- map_trait(comm, x, FUN = sd, shp=poly)

par(mar=rep(0,4))
plot_swatch(y, y$traits, k=20, border=NA, leg=3)

## -----------------------------------------------------------------------------
data(africa)
sparse_comm <- africa$comm

tree <- africa$phylo
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(sparse_comm)))
pb <- phylobeta(sparse_comm, tree)

## ---- message=FALSE, results='hide', warning=FALSE----------------------------
y <- phyloregion(pb[[1]], shp=africa$polys)

## -----------------------------------------------------------------------------
plot_NMDS(y, cex=3)
text_NMDS(y)

par(mar=rep(0,4))
plot(y, palette="NMDS")

## ---- eval=TRUE---------------------------------------------------------------
sessionInfo()

