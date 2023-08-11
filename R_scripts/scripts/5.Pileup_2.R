# Benjamin Jean-Marie Tremblay

library(Biostrings)
library(ape)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)

# Missing from plasmid:
#
# Punta_Candelero−Tooth
# Peru−NA40−Bone
# Abusir1607−Tooth
# 407.86
# 1586−Z1
# ATCC−19406

gffF <- "data/E88_plasmid.gff3"
gff0 <- import(gffF)
gff <- suppressMessages(readr::read_tsv(gffF, col_names = FALSE, comment = "#"))
gff2 <- gff[gff[[2]] == "RefSeq", ]
gff2 <- gff2[, c(1, 4:5, 3)]
colnames(gff2) <- c("chr", "start", "end", "feature")
gff2 <- gff2[!is.na(gff2$chr), ]
genes <- gff2[gff2$feature == "gene", ]
genesOld <- genes
genes <- as.data.frame(reduce(GRanges(genes)))

gffCF <- "data/E88_chromosome.gff3"
gffC0 <- import(gffCF)
gffC <- suppressMessages(readr::read_tsv(gffCF, col_names = FALSE, comment = "#"))
gffC2 <- gffC[gffC[[2]] == "RefSeq", ]
gffC2 <- gffC2[, c(1, 4:5, 3)]
colnames(gffC2) <- c("chr", "start", "end", "feature")
gffC2 <- gffC2[!is.na(gffC2$chr), ]
genesC <- gffC2[gffC2$feature == "gene", ]
genesCOld <- genesC
genesC <- as.data.frame(reduce(GRanges(genesC)))

# chr: 1-2799251
# plasmid: 2799252-2873333

msa <- readDNAStringSet("data/clean.full.reduced.aln")
tree <- read.tree("data/final-snippy-fasttree-oct23.newick")
tree$tip.label <- gsub("'", "", tree$tip.label, fixed = TRUE)

msaM <- strsplit(as.character(msa), "", TRUE)
msaM <- as.matrix(t(list2DF(msaM)))

ref <- msaM["Reference", ]

# gap: NA
# same as ref: Ref
# snp: Alt

msaMsnp <- msaM

for (i in 1:nrow(msaM)) {

  msaMsnp[i, msaMsnp[i, ] == ref] <- "Ref"
  msaMsnp[i, msaMsnp[i, ] == "-"] <- NA
  msaMsnp[i, !is.na(msaMsnp[i, ]) & msaMsnp[i, ] != "Ref"] <- "Alt"

  msaM[i, msaM[i, ] == "-"] <- NA
  msaM[i, !is.na(msaM[i, ])] <- "Ref"

}

msaMsnpchr <- msaMsnp[, 1:2799251]
msaMsnppld <- msaMsnp[, 2799252:2873333]

msaMchr <- msaM[, 1:2799251]
msaMpld <- msaM[, 2799252:2873333]

pctId <- numeric(nrow(msaM))
pctIdGapless <- numeric(nrow(msaM))
pctIdChr <- numeric(nrow(msaM))
pctIdPld <- numeric(nrow(msaM))
pctIdChrGapless <- numeric(nrow(msaM))
pctIdPldGapless <- numeric(nrow(msaM))

for (i in 1:nrow(msaM)) {
  pctId[i] <- sum(!is.na(msaMsnp[i, ]) & msaMsnp[i, ] == "Ref") / ncol(msaMsnp)
  pctIdGapless[i] <- sum(!is.na(msaMsnp[i, ]) & msaMsnp[i, ] == "Ref") / sum(!is.na(msaMsnp[i, ]))
  pctIdChr[i] <- sum(!is.na(msaMsnpchr[i, ]) & msaMsnpchr[i, ] == "Ref") / ncol(msaMsnpchr)
  pctIdPld[i] <- sum(!is.na(msaMsnppld[i, ]) & msaMsnppld[i, ] == "Ref") / ncol(msaMsnppld)
  pctIdChrGapless[i] <- sum(!is.na(msaMsnpchr[i, ]) & msaMsnpchr[i, ] == "Ref") / sum(!is.na(msaMsnpchr[i, ]))
  pctIdPldGapless[i] <- sum(!is.na(msaMsnppld[i, ]) & msaMsnppld[i, ] == "Ref") / sum(!is.na(msaMsnppld[i, ]))
}

names(pctId) <- rownames(msaM)
names(pctIdGapless) <- rownames(msaM)
names(pctIdChr) <- rownames(msaM)
names(pctIdPld) <- rownames(msaM)
names(pctIdChrGapless) <- rownames(msaM)
names(pctIdPldGapless) <- rownames(msaM)

pctId <- pctId * 100
pctIdGapless <- pctIdGapless * 100
pctIdChr <- pctIdChr * 100
pctIdPld <- pctIdPld * 100
pctIdChrGapless <- pctIdChrGapless * 100
pctIdPldGapless <- pctIdPldGapless * 100

pctId[is.na(pctId)] <- 0
pctIdGapless[is.na(pctIdGapless)] <- 0
pctIdChr[is.na(pctIdChr)] <- 0
pctIdPld[is.na(pctIdPld)] <- 0
pctIdChrGapless[is.na(pctIdChrGapless)] <- 0
pctIdPldGapless[is.na(pctIdPldGapless)] <- 0

pld_tent <- c(68640, 72587)
msaTent <- msaMsnp[, 2799252:2873333][, pld_tent[1]:pld_tent[2]]
pctIdTent <- apply(msaTent, 1, function(x) {x=x[!is.na(x)];sum(x=="Ref")/length(x)})



makeGRange <- function(x, name) {
  x <- Rle(x)
  GRanges(name, IRanges(cumsum(c(0, runLength(x)[-nrun(x)])),
      width = runLength(x)), base = runValue(x))
}

chrGR <- vector("list", nrow(msaM))
pldGR <- vector("list", nrow(msaM))

for (i in 1:nrow(msaM)) {
  chrGR[[i]] <- cbind(as.data.frame(makeGRange(msaMchr[i, ], "chr")), sample = rownames(msaM)[i])
  pldGR[[i]] <- cbind(as.data.frame(makeGRange(msaMpld[i, ], "plasmid")), sample = rownames(msaM)[i])
}

chrDf <- do.call(rbind, chrGR)[, c("start", "end", "base", "sample")]
pldDf <- do.call(rbind, pldGR)[, c("start", "end", "base", "sample")]

chrDf <- chrDf[!is.na(chrDf$base), ]
pldDf <- pldDf[!is.na(pldDf$base), ]

chrDf$PercentID <- pctIdChrGapless[chrDf$sample]
pldDf$PercentID <- pctIdPldGapless[pldDf$sample]

genesPld <- data.frame(row.names = NULL, start = genes$start, end = genes$end,
  base = "Ref", sample = "Genes", PercentID = NA_real_)
genesChr <- data.frame(row.names = NULL, start = genesC$start, end = genesC$end,
  base = "Ref", sample = "Genes", PercentID = NA_real_)

TeNT <- data.frame(row.names = NULL, start = 68640, end = 72587, base = "Ref",
  sample = factor("Genes", levels = c("Genes", tree$tip.label)))
colT <- data.frame(row.names = NULL, start = 39438, end = 42413, base = "Ref",
  sample = factor("Genes", levels = c("Genes", tree$tip.label)))
repA <- data.frame(row.names = NULL, start = 20250, end = 21941, base = "Ref",
  sample = factor("Genes", levels = c("Genes", tree$tip.label)))

chrDf2 <- rbind(genesChr, chrDf)
pldDf2 <- rbind(genesPld, pldDf)

chrDf2$sample <- factor(chrDf2$sample, levels = c("Genes", tree$tip.label))
pldDf2$sample <- factor(pldDf2$sample, levels = c("Genes", tree$tip.label))

p1 <- ggplot(pldDf2, aes(x = start, xend = end, y = 0, yend = 0, colour = PercentID)) +
  geom_segment(lwd = 2.2, lineend = "butt") +
  scale_colour_viridis_c(na.value = "black") +
  geom_segment(data = TeNT, aes(x = start, xend = end, y = 0, yend = 0),
    colour = "red", lwd = 2.2, lineend = "butt") +
  geom_segment(data = colT, aes(x = start, xend = end, y = 0, yend = 0),
    colour = "blue", lwd = 2.2, lineend = "butt") +
  geom_segment(data = repA, aes(x = start, xend = end, y = 0, yend = 0),
    colour = "purple", lwd = 2.2, lineend = "butt") +
  facet_wrap(~sample, ncol = 1, strip.position = "right", drop = FALSE) +
  theme_void() +
  theme(
    panel.spacing = unit(0.1, "mm"),
    strip.text.y.right = element_text(hjust = 0, size = unit(6, "pt"))
  )

p2 <- ggplot(chrDf2, aes(x = start, xend = end, y = 0, yend = 0, colour = PercentID)) +
  geom_segment(lwd = 2.2, lineend = "butt") +
  scale_colour_viridis_c(na.value = "black") +
  facet_wrap(~sample, ncol = 1, strip.position = "right") +
  theme_void() +
  theme(
    panel.spacing = unit(0.1, "mm"),
    strip.text.y.right = element_text(hjust = 0, size = unit(6, "pt"))
  )

ggsave("figures/Plasmid_pctId_oct24.pdf", p1, width = 9.8, height = 6.12)
ggsave("figures/Chromosome_pctId_oct23.pdf", p2, width = 9.8, height = 6.12)

pldGR <- makeGRangesFromDataFrame(pldDf2, seqnames.field = "base", keep = TRUE)

TeNTGR <- pldGR[overlapsAny(pldGR, GRanges("Ref", IRanges(68640, 72587)))]
start(TeNTGR)[start(TeNTGR) < 68640] <- 68640
end(TeNTGR)[end(TeNTGR) > 72587] <- 72587

colTGR <- pldGR[overlapsAny(pldGR, GRanges("Ref", IRanges(39438, 42413)))]
start(colTGR)[start(colTGR) < 39438] <- 39438
end(colTGR)[end(colTGR) > 42413] <- 42413

repAGR <- pldGR[overlapsAny(pldGR, GRanges("Ref", IRanges(20250, 21941)))]
start(repAGR)[start(repAGR) < 20250] <- 20250
end(repAGR)[end(repAGR) > 21941] <- 21941

p3 <- ggplot(pldDf2, aes(x = start, xend = end, y = 0, yend = 0)) +
  geom_segment(lwd = 2.2, lineend = "butt", colour = "black") +
  geom_segment(data = as.data.frame(TeNTGR),
    aes(x = start, xend = end, y = 0, yend = 0),
    colour = "red", lwd = 2.2, lineend = "butt") +
  geom_segment(data = as.data.frame(colTGR),
    aes(x = start, xend = end, y = 0, yend = 0),
    colour = "blue", lwd = 2.2, lineend = "butt") +
  geom_segment(data = as.data.frame(repAGR),
    aes(x = start, xend = end, y = 0, yend = 0),
    colour = "purple", lwd = 2.2, lineend = "butt") +
  facet_wrap(~sample, ncol = 1, strip.position = "right", drop = FALSE) +
  theme_void() +
  theme(
    panel.spacing = unit(0.1, "mm"),
    strip.text.y.right = element_text(hjust = 0, size = unit(6, "pt"))
  )

pldDf3 <- pldDf2[pldDf2$sample != "Genes", ]
pldDf3$PercentID <- 0
pldDf3$sam2 <- as.character(pldDf3$sample)

chrDf3 <- chrDf2[chrDf2$sample != "Genes", ]
chrDf3$PercentID <- 0
chrDf3$sam2 <- as.character(chrDf3$sample)

calc_pctid <- function(x, x1, x2, s) {
  x <- x[s, x1:x2]
  100 * (sum(!is.na(x) & x == "Ref") / sum(!is.na(x)))
}

for (i in 1:nrow(chrDf3)) {
  if (i %% 1000 == 0) message(i, " / ", nrow(chrDf3))
  chrDf3$PercentID[i] <- calc_pctid(msaMsnpchr, chrDf3$start[i]+1, chrDf3$end[i], chrDf3$sam2[i])
}
for (i in 1:nrow(pldDf3)) {
  if (i %% 1000 == 0) message(i, " / ", nrow(pldDf3))
  pldDf3$PercentID[i] <- calc_pctid(msaMsnppld, pldDf3$start[i]+1, pldDf3$end[i], pldDf3$sam2[i])
}

p4 <- ggplot(chrDf3, aes(x = start, xend = end, y = 0, yend = 0,colour=PercentID)) +
  geom_segment(lwd = 2.2, lineend = "butt")+
  facet_wrap(~sample, ncol = 1, strip.position = "right") +
  scale_colour_viridis_c(na.value="grey50") +
  theme_void() +
  theme(
    panel.spacing = unit(0.1, "mm"),
    strip.text.y.right = element_text(hjust = 0, size = unit(6, "pt"))
  )

ggsave("Plasmid_simple_nov7.pdf", p3, width = 9.12, height = 6.12)
ggsave("Chromosome_pctid_nov27.pdf", p4, width = 9.12, height = 6.12)



###############################################################################

msa <- readDNAStringSet("/data/consensus_variants.E88_aligned.noInsertions.fastANDREW-a")

translate2 <- function(x) {
  cdns <- universalmotif::window_string(as.character(x), 3, 0)
  aa <- character(length(cdns))
  for (i in seq_along(aa)) {
    aa[i] <- tryCatch(as.character(translate(DNAString(cdns[i]))), error = function(e) "-")
  }
  AAString(paste0(aa, collapse = ""))
}
translate3 <- function(x) {
  y <- vector("list", length(x))
  for (i in seq_along(y)) {
    message(i, " / ", length(y))
    y[[i]] <- translate2(x[i])
  }
  unname(y)
}

msa_aa <- translate3(msa)
msa_aa <- AAStringSet(msa_aa)
names(msa_aa) <- names(msa)
writeXStringSet(msa_aa, "tent_seqs.faa")

calc_pct <- function(x, y) {
  x <- strsplit(as.character(x), "", TRUE)[[1]]
  y <- strsplit(as.character(y), "", TRUE)[[1]]
  gaps <- x == "-" | y == "-"
  x <- x[!gaps]
  y <- y[!gaps]
  100 * (sum(x == y) / length(x))
}

calc_cov <- function(x) {
  x <- strsplit(as.character(x), "", TRUE)[[1]]
  100 * (sum(x != "-") / length(x))
}

isRef <- !grepl("^SAM", names(msa))

id_all_v_all <- matrix(0, nrow = sum(isRef), ncol = sum(!isRef),
  dimnames = list(names(msa)[isRef], names(msa)[!isRef]))
for (i in 1:ncol(id_all_v_all)) {
  for (j in 1:nrow(id_all_v_all)) {
    id_all_v_all[j, i] <- calc_pct(msa_aa[[rownames(id_all_v_all)[j]]],
      msa_aa[[colnames(id_all_v_all)[i]]])
  }
}

id_all_v_all <- reshape2::melt(id_all_v_all)
id_all_v_all[[1]] <- as.character(id_all_v_all[[1]])
id_all_v_all[[2]] <- as.character(id_all_v_all[[2]])
colnames(id_all_v_all) <- c("ModernTeNT", "acTeNT", "PctID")

ids <- data.frame(row.names = NULL,
  BioSample = names(msa)[grepl("^SAM", names(msa))],
  Aln_cov = NA_real_,
  PctID_E88 = NA_real_,
  PctID_highest = NA_real_,
  Which_highest = NA_character_
)

for (i in 1:nrow(ids)) {
  ids$Aln_cov[i] <- calc_cov(msa_aa[[ids$BioSample[i]]])
  ids$PctID_E88[i] <- calc_pct(msa_aa[["E88_Reference"]], msa_aa[[ids$BioSample[i]]])
  ids$PctID_highest[i] <- max(id_all_v_all$PctID[id_all_v_all$acTeNT==ids$BioSample[i]])
  ids$Which_highest[i] <- paste0(
    id_all_v_all$ModernTeNT[id_all_v_all$acTeNT==ids$BioSample[i]][
      id_all_v_all$PctID[id_all_v_all$acTeNT==ids$BioSample[i]] == 
      max(id_all_v_all$PctID[id_all_v_all$acTeNT==ids$BioSample[i]])
    ], collapse = ", ")
}


id_mod_v_mod <- matrix(0, nrow = sum(isRef), ncol = sum(isRef),
  dimnames = list(names(msa)[isRef], names(msa)[isRef]))
for (i in 1:ncol(id_mod_v_mod)) {
  for (j in 1:nrow(id_mod_v_mod)) {
    id_mod_v_mod[j, i] <- calc_pct(msa_aa[[rownames(id_mod_v_mod)[j]]],
      msa_aa[[colnames(id_mod_v_mod)[i]]])
  }
}

id_mod_v_mod <- reshape2::melt(id_mod_v_mod)
id_mod_v_mod[[1]] <- as.character(id_mod_v_mod[[1]])
id_mod_v_mod[[2]] <- as.character(id_mod_v_mod[[2]])
colnames(id_mod_v_mod) <- c("ModernTeNT", "ModernTeNT", "PctID")

readr::write_tsv(ids, "results/TeNT_PctID_supp.txt")
readr::write_tsv(id_all_v_all, "results/TeNT_PctID_ancient.txt")
readr::write_tsv(id_mod_v_mod, "results/TeNT_PctID_modern.txt")

