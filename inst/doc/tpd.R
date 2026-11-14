## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4, message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(phyloregion)
library(terra)

data(africa)
p <- vect(system.file("ex/sa.json", package = "phyloregion"))

tr <- data.frame(species = africa$trait$species,
                 height = log10(africa$trait$h_max),
                 habit = africa$trait$habit)
head(tr)

## -----------------------------------------------------------------------------
fa <- tpd_alpha(africa$comm, tr)
head(fa)

## ----fig.height = 3.2---------------------------------------------------------
m <- merge(p, fa, by = "grids")
par(mfrow = c(1, 3), mar = c(1, 1, 2, 3))
plot(m, "FRic", main = "FRic", axes = FALSE, border = NA,
     col = hcl.colors(20, "YlOrRd", rev = TRUE))
plot(m, "FEve", main = "FEve", axes = FALSE, border = NA,
     col = hcl.colors(20, "YlOrRd", rev = TRUE))
plot(m, "FRed", main = "FRed", axes = FALSE, border = NA,
     col = hcl.colors(20, "YlOrRd", rev = TRUE))

## -----------------------------------------------------------------------------
fb <- tpd_beta(africa$comm, tr)
sapply(fb, function(d) round(mean(d), 3))

## ----fig.height = 4.5---------------------------------------------------------
keep <- fa$grids[fa$SR >= 10]
fbk <- as.dist(as.matrix(fb$beta.total)[keep, keep])
fr <- phyloregion(fbk, k = 5, pol = p)
plot(fr)

## -----------------------------------------------------------------------------
status <- ifelse(africa$IUCN$IUCN %in% c("CR", "EN", "VU"),
                 "threatened", "not_threatened")
names(status) <- africa$IUCN$Species
g <- tpd_groups(status)
rowSums(g)

tpd_alpha(g, tr)[, c("grids", "SR", "FRic", "FRed", "CWM_height")]
round(as.matrix(tpd_beta(g, tr)$beta.total), 3)

## -----------------------------------------------------------------------------
pr <- c(LC = 0.001, NT = 0.01, VU = 0.1, EN = 0.5, CR = 0.9)
prob <- setNames(pr[africa$IUCN$IUCN], africa$IUCN$Species)
set.seed(1)
ex <- tpd_extinction(tr, prob, nsim = 50)
ex$summary[ex$summary$quantity %in% c("lost", "FRic", "overlap"), ]

## -----------------------------------------------------------------------------
sessionInfo()

