# Benjamin Jean-Marie Tremblay
# 2021-07-09

library(Biostrings)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)

msa <- readDNAStringSet("data/core.full.cleaned.aln")
msa <- strsplit(as.character(msa), "", TRUE)
msa <- as.matrix(t(list2DF(msa)))

msaNA <- msa
msaNA[msaNA != "-"] <- "B"
msaNA[msaNA != "B"] <- NA

msaRef <- msa["Reference", , drop = FALSE]
msaRef[msaRef == "-"] <- NA

msaNoRef <- msa[rownames(msa) != "Reference", ]
for (i in 1:nrow(msaNoRef)) msaNoRef[i, ][msaNoRef[i, ] == msaRef[1, ]] <- "Ref"
msaNoRef[msaNoRef != "Ref" & msaNoRef != "-"] <- "Alt"
msaNoRef[msaNoRef == "-"] <- NA

msaRef[!is.na(msaRef)] <- "Ref"

makeGRange <- function(x, name) {
  x <- Rle(x)
  GRanges(name, IRanges(cumsum(c(0, runLength(x)[-nrun(x)])),
      width = runLength(x)), base = runValue(x))
}

msaNAGR <- structure(vector("list", nrow(msaNA)), names = names(msaNA))
for (i in 1:nrow(msaNA)) {
  msaNAGR[[i]] <- as.data.frame(makeGRange(msaNA[i, ], "theseq"))
  msaNAGR[[i]] <- msaNAGR[[i]][!is.na(msaNAGR[[i]]$base), ]
}
names(msaNAGR) <- rownames(msaNA)

plasmidNoRef <- msaNoRef[, 2799251:ncol(msaNoRef)]
plasmidRef <- msaRef[, 2799251:ncol(msaRef), drop = FALSE]

plasmidNAGR <- lapply(msaNAGR, function(x) x[x$start > 2799250, ])

plasmidNAGR <- lapply(plasmidNAGR, function(x) {
  x$start <- x$start - 2799251
  x$end <- x$end - 2799251
  x
})

chromNA <- msaNA[, 1:2799250]
chromNoRef <- msaNoRef[, 1:2799250]
chromRef <- msaRef[, 1:2799250, drop = FALSE]
chromNAGR <- lapply(msaNAGR, function(x) x[x$start < 2799251, ])

fixEmpty <- function(x) {
  if (!nrow(x))
    data.frame(seqnames = "theseq", start = 0, end = 1, width = 1, base = "E")
  else x
}

chromNAGR <- lapply(chromNAGR, fixEmpty)
plasmidNAGR <- lapply(plasmidNAGR, fixEmpty)

msaNoRefPctID <- structure(numeric(nrow(msaNoRef)), names = rownames(msaNoRef))
chromNoRefPctID <- structure(numeric(nrow(chromNoRef)), names = rownames(chromNoRef))
plasmidNoRefPctID <- structure(numeric(nrow(plasmidNoRef)), names = rownames(plasmidNoRef))

calcPctID <- function(x) {
  x <- x[!is.na(x)]
  (sum(x == "Ref") / length(x)) * 100
}
for (i in seq_along(msaNoRefPctID))
  msaNoRefPctID[i] <- calcPctID(msaNoRef[i, ])
for (i in seq_along(chromNoRefPctID))
  chromNoRefPctID[i] <- calcPctID(chromNoRef[i, ])
for (i in seq_along(plasmidNoRefPctID))
  plasmidNoRefPctID[i] <- calcPctID(plasmidNoRef[i, ])

CircosOrder <- data.frame(
  BioSample = names(msaNoRefPctID)[grepl("SAM", names(msaNoRefPctID))],
  PercentID = unname(msaNoRefPctID)[grepl("SAM", names(msaNoRefPctID))],
  ChromPercentID = unname(chromNoRefPctID)[grepl("SAM", names(chromNoRefPctID))],
  PlasmidPercentID = unname(plasmidNoRefPctID)[grepl("SAM", names(plasmidNoRefPctID))]
)
CircosOrder <- CircosOrder[order(CircosOrder$PercentID, decreasing = TRUE), ]
readr::write_tsv(CircosOrder, "results/BioSampleCircosOrder.tsv")

plasmidNoRefPctID <- plasmidNoRefPctID[rev(CircosOrder$BioSample)]
chromNoRefPctID <- chromNoRefPctID[rev(CircosOrder$BioSample)]

plasmidSAMPctID <- plasmidNoRefPctID[grepl("SAM", names(plasmidNoRefPctID))]
chromSAMPctID <- chromNoRefPctID[grepl("SAM", names(chromNoRefPctID))]

plasmidSAMGR <- plasmidNAGR[names(plasmidSAMPctID)]
chromSAMGR <- chromNAGR[names(chromSAMPctID)]

plasmidSAMPctIDcut <- as.integer(cut(plasmidSAMPctID, 89:100, include.lowest = TRUE))
chromSAMPctIDcut <- as.integer(cut(chromSAMPctID, 89:100, include.lowest = TRUE))

PctIDcolours <- rev(viridis::viridis_pal()(12))
PctIDlegend <- Legend(at = rev(1:10), legend_gp = gpar(fill = rev(PctIDcolours)),
  labels = rev(c("89-90%", "90-91%", "91-92%", "92-93%", "93-94%", "94-95%",
      "95-96%", "96-97%", "97-98%", "98-99%", "99-100%")),
  grid_height = unit(4, "mm"), grid_width = unit(4, "mm"),
  gap = unit(0, "mm"), labels_gp = gpar(fontsize = 10))

ideogramChrom <- data.frame("theseq", 0, 2799250, "A", "gneg")
ideogramPlasmid <- data.frame("theseq", 0, 74081, "A", "gneg")

genesPlasmid <- readr::read_tsv("data/E88_plasmid.gff3", comment = "#",
  col_names = FALSE)
genesChrom <- readr::read_tsv("data/E88_chromosome.gff3", comment = "#",
  col_names = FALSE)

genesPlasmid <- genesPlasmid[genesPlasmid$X2 == "RefSeq" & genesPlasmid$X3 == "gene", ]
genesPlasmid$X5[nrow(genesPlasmid)] <- 74081
genesPlasmidDf <- data.frame(
  seqnames = "theseq", start = genesPlasmid$X4, end = genesPlasmid$X5,
  width = genesPlasmid$X5 - genesPlasmid$X4 + 1, strand = "*"
)

genesChrom <- genesChrom[genesChrom$X2 == "RefSeq" & genesChrom$X3 == "gene", ]
genesChromDf <- data.frame(
  seqnames = "theseq", start = genesChrom$X4, end = genesChrom$X5,
  width = genesChrom$X5 - genesChrom$X4 + 1, strand = "*"
)

circos.clear()
pdf("figures/CircosPercentID_chromosome_SAM.pdf")

circos.par(
  "start.degree" = 45,
  cell.padding = c(0, 0, 0, 0),
  gap.degree = 45.5,
  points.overflow.warning = FALSE,
  track.margin = c(0.0025, 0)
)

circos.initializeWithIdeogram(ideogramChrom,
  plotType = NULL,
  axis.labels.cex = 0.8,
  major.by = 500000
)

circos.track(ylim = c(0, 1), track.height = 0.001,
  panel.fun = function(x, y) {
    circos.genomicAxis(
      major.at = c(.5e6, 1e6, 1.5e6, 2e6, 2.5e6),
      labels.cex = 0.9,
      labels = c("0.5MB", "1.0MB", "1.5MB", "2.0MB", "2.5MB"))
})

circos.genomicTrackPlotRegion(genesChromDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.01,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "white", border = NA)
})
circos.genomicTrackPlotRegion(genesChromDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.02,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "black", border = NA)
})
circos.genomicTrackPlotRegion(genesChromDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.01,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "white", border = NA)
})

for (i in seq_along(chromSAMGR)) {
  message(i, "/", length(chromSAMGR))
  circos.genomicTrackPlotRegion(chromSAMGR[[i]],
    ylim = c(0, 1), bg.border = NA, track.height = 0.02,
    panel.fun = function(region, value, ...) {
      if (i == 4) {
        circos.text(2799250 + convert_x(3, "mm"), 0.4, 35, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 9) {
        circos.text(2799250 + convert_x(3, "mm"), 0.4, 30, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 14) {
        circos.text(2799250 + convert_x(3, "mm"), 0.4, 25, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 19) {
        circos.text(2799250 + convert_x(3, "mm"), 0.4, 20, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 24) {
        circos.text(2799250 + convert_x(2.5, "mm"), 0.5, 15, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 29) {
        circos.text(2799250 + convert_x(2.5, "mm"), 0.5, 10, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 34) {
        circos.text(2799250 + convert_x(3, "mm"), 0.6, 5, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      circos.text(2799250 + 3779, 0.7, "-", cex = 0.8,
        facing = "downward", adj = c(0, 0.5))
      circos.genomicRect(region, ytop = 1, ybottom = 0,
        col = PctIDcolours[chromSAMPctIDcut[i]], border = NA)
    })
}

draw(PctIDlegend, x = unit(0.06, "npc"), y = unit(0.15, "npc"))

dev.off()

circos.clear()
pdf("figures/CircosPercentID_plasmid_SAM.pdf")

circos.par(
  "start.degree" = 45,
  cell.padding = c(0, 0, 0, 0),
  gap.degree = 45.5,
  points.overflow.warning = FALSE,
  track.margin = c(0.0025, 0)
)

circos.initializeWithIdeogram(ideogramPlasmid,
  plotType = NULL,
  axis.labels.cex = 0.8,
  major.by = 10000
)
circos.track(ylim = c(0, 1), track.height = 0.001,
  panel.fun = function(x, y) {
    circos.genomicAxis(
      major.at = c(1e4, 3e4, 5e4, 7e4),
      labels = c("10KB", "30KB", "50KB", "70KB"),
      labels.cex = 0.9)
})

circos.genomicTrackPlotRegion(genesPlasmidDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.01,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "white", border = NA)
})
circos.genomicTrackPlotRegion(genesPlasmidDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.02,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "black", border = NA)
    circos.genomicRect(genesPlasmidDf[83, ], ytop = 1, ybottom = 0, col = "red", border = NA)
})
circos.genomicTrackPlotRegion(genesPlasmidDf,
  ylim = c(0, 1), bg.border = NA, track.height = 0.01,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ytop = 1, ybottom = 0, col = "white", border = NA)
})

for (i in seq_along(plasmidSAMGR)) {
  message(i, "/", length(plasmidSAMGR))
  circos.genomicTrackPlotRegion(plasmidSAMGR[[i]],
    ylim = c(0, 1), bg.border = NA, track.height = 0.02,
    panel.fun = function(region, value, ...) {
      if (i == 4) {
        circos.text(74081 + convert_x(3, "mm"), 0.4, 35, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 9) {
        circos.text(74081 + convert_x(3, "mm"), 0.4, 30, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 14) {
        circos.text(74081 + convert_x(3, "mm"), 0.4, 25, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 19) {
        circos.text(74081 + convert_x(3, "mm"), 0.4, 20, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 24) {
        circos.text(74081 + convert_x(2.5, "mm"), 0.5, 15, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 29) {
        circos.text(74081 + convert_x(2.5, "mm"), 0.5, 10, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      if (i == 34) {
        circos.text(74081 + convert_x(3, "mm"), 0.6, 5, cex = 0.9,
          facing = "downward", adj = c(0, 0.5))
      }
      circos.text(74081 + 100, 0.7, "-", cex = 0.8,
        facing = "downward", adj = c(0, 0.5))
      circos.genomicRect(region, ytop = 1, ybottom = 0,
        col = PctIDcolours[plasmidSAMPctIDcut[i]], border = NA)
    })
}

draw(PctIDlegend, x = unit(0.06, "npc"), y = unit(0.15, "npc"))

dev.off()
