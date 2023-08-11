# Benjamin Jean-Marie Tremblay
# 2021-07-10

library(ggplot2)

tree <- ape::read.tree(text=readr::read_lines("data/tent-tree-new")[1])

msa <- Biostrings::readDNAStringSet("data/consensus_variants.E88_aligned.noInsertions.fastANDREW-a")
names(msa)[names(msa) == "E88_Reference"] <- "Reference"
msa <- as.character(msa)
msa <- t(sapply(msa, \(.) strsplit(., "", TRUE)[[1]]))

rownames(msa)[rownames(msa) == "1586-U1"] <- "U1"

Ref <- msa["Reference", ]
msaOriginal <- msa

for (i in 1:ncol(msa)) {
  msa[, i][msa[, i] == Ref[i]] <- "Black"
  msa[, i][!msa[, i] %in% c("Black", "-")] <- "Red"
  msa[, i][msa[, i] == "-"] <- NA
}

SAMs <- grepl("SAM", rownames(msa))

for (i in 1:ncol(msa)) {
  Ancient <- msa[SAMs, i]
  New <- msa[!SAMs, i]
  New <- New[!is.na(New)]
  if (any(New == "Red")) next
  if (any(Ancient[!is.na(Ancient)] == "Red")) {
    Ancient[Ancient == "Red"] <- "Yellow"
    msa[SAMs, i] <- Ancient
  }
}

rowOrder <- readr::read_lines("data/ToxinOrder.txt")
rowOrder <- rownames(msa)[pmatch(rowOrder, rownames(msa))]
msa <- msa[rowOrder, ]

SNPcount <- apply(msa, 1, \(x) sum(x == "Yellow", na.rm = TRUE))
SNPcount <- data.frame(
  BioSample = names(SNPcount),
  AncientSNPsCount = unname(SNPcount)
)
readr::write_tsv(SNPcount, "results/ToxinAncientSNPsCount.tsv")

msaDF <- reshape2::melt(msa)
msaDF$Var1 <- factor(msaDF$Var1, levels = rev(rownames(msa)))
msaDF$Var2 <- as.integer(msaDF$Var2)
msaDF <- msaDF[!is.na(msaDF$value), ]

msaDFsnps <- msaDF[msaDF$value != "Black", ]
colnames(msaDFsnps) <- c("Sample", "Position", "IsAncient")
msaDFsnps$IsAncient[msaDFsnps$IsAncient == "Red"] <- "No"
msaDFsnps$IsAncient[msaDFsnps$IsAncient == "Yellow"] <- "Yes"
msaDFsnps$Ref <- Ref[msaDFsnps$Position]
msaDFsnps$Alt <- NA_character_
msaDFsnps$Sample <- as.character(msaDFsnps$Sample)
for (i in 1:nrow(msaDFsnps)) {
  msaDFsnps$Alt[i] <- msaOriginal[msaDFsnps$Sample[i], msaDFsnps$Position[i]]
}

readr::write_tsv(msaDFsnps, "results/ToxinSNPs.tsv")
msaDFsnps2 <- dplyr::filter(msaDFsnps, IsAncient == "Yes" & grepl("^SAM", Sample))
msaDFsnps2$Sample <- gsub(".consensus.masked", "", msaDFsnps2$Sample, fixed = TRUE)
msaDFsnps2 <- msaDFsnps2[, -3]
msaDFsnps2 <- aggregate(msaDFsnps2[1], msaDFsnps2[-1],
  FUN = \(x) paste(sort(unique(x)), collapse = ", "))
colnames(msaDFsnps2) <- c("Ancient SNP", "Ref", "Alt", "Samples")
msaDFsnps2 <- msaDFsnps2[order(msaDFsnps2[[1]]), ]
readr::write_tsv(msaDFsnps2, "results/ToxinAncientSNPsUnique.tsv")

msaDF_tmp <- msaDF[msaDF$value %in% c("Red", "Yellow"), ]
msaDF_tmp$Var2 <- msaDF_tmp$Var2 + 1
msaDF <- rbind(msaDF_tmp, msaDF)
msaDF_tmp$Var2 <- msaDF_tmp$Var2 + 1
msaDF <- rbind(msaDF_tmp, msaDF)
msaDF_tmp$Var2 <- msaDF_tmp$Var2 - 3
msaDF <- rbind(msaDF_tmp, msaDF)
msaDF_tmp$Var2 <- msaDF_tmp$Var2 - 1
msaDF <- rbind(msaDF_tmp, msaDF)

msaDF <- msaDF[order(match(msaDF$value, c("Black", "Red", "Yellow"))), ]

ordNew <- gsub("'", "", tree$tip.label, fixed=T)
ordNew[c(1, 3, 5)] <- paste0(ordNew[c(1, 3, 5)], "d")
msaDF_tmp <- msaDF[as.character(msaDF[[1]]) %in% ordNew, ]

p <- ggplot(msaDF_tmp, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0),
    breaks = c(1000, 2000, 3000),
    labels = c("1KB", "2KB", "3KB")) +
  scale_y_discrete(position = "right", limits = rev(ordNew)) +
  scale_fill_manual(
    values = c(Black = "black", Yellow = "yellow", Red = "red"),
    na.value = "white") +
  ylab(element_blank()) +
  xlab(element_blank()) +
  theme_bw() +
  theme(panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), unit = "cm"),
    axis.text.y = element_text(colour = "black", size = rel(0.5)),
    axis.text.x = element_text(colour = "black"))
ggsave("figures/ToxinSNPs2.png",
 plot = p, width = 6.75, height = 4.2)
ggsave("figures/ToxinSNPs2.pdf",
 plot = p, width = 6.75, height = 4.2)



###############################################################################
# Create a lollipop chart of the ancient SNPs (so using msaDFsnps)

library(S4Vectors)

msaDF3 <- reshape2::melt(msa)
msaDF3$Var1 <- factor(msaDF3$Var1, levels = rev(rownames(msa)))
msaDF3$Var2 <- as.integer(msaDF3$Var2)
msaDF3 <- msaDF3[!is.na(msaDF3$value), ]

msaDFsnps3 <- msaDF3
colnames(msaDFsnps3) <- c("Sample", "Position", "IsAncient")
msaDFsnps3$IsAncient[msaDFsnps3$IsAncient == "Red"] <- "No"
msaDFsnps3$IsAncient[msaDFsnps3$IsAncient == "Yellow"] <- "Yes"
msaDFsnps3$IsAncient[msaDFsnps3$IsAncient == "Black"] <- NA
msaDFsnps3$Ref <- Ref[msaDFsnps3$Position]
msaDFsnps3$Alt <- NA_character_
msaDFsnps3$Sample <- as.character(msaDFsnps3$Sample)
for (i in 1:nrow(msaDFsnps3)) {
  msaDFsnps3$Alt[i] <- msaOriginal[msaDFsnps3$Sample[i], msaDFsnps3$Position[i]]
}

msaDFsnps3$y <- ifelse(!is.na(msaDFsnps3$IsAncient) & msaDFsnps3$IsAncient == "Yes", 1, 0)
msaDFsnps3$lab <- paste0(msaDFsnps3$Ref, msaDFsnps3$Position, msaDFsnps3$Alt)
msaDFsnps3$lab[is.na(msaDFsnps3$IsAncient) | msaDFsnps3$IsAncient != "Yes"] <- NA

msaDFsnps3$Sample <- gsub(".consensus.masked", "", msaDFsnps3$Sample, fixed = TRUE)

msaDFsnps3 <- msaDFsnps3[msaDFsnps3$Sample %in%
  unique(msaDFsnps3$Sample[msaDFsnps3$y == 1]), ]

msaDFsnpsL <- tapply(as.integer(msaDFsnps3$y == 0) + 1,
  msaDFsnps3$Sample, function(x) list(runLength(Rle(x))))

msaDFsnps3$PosM1 <- msaDFsnps3$Position - 1

ggplot(msaDFsnps3[msaDFsnps3$y > 0, ],
  aes(x = Position, y = y)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = y), lwd = 0.5) +
  geom_point(shape = 21, fill = "white", colour = "black", size = 2) +
  geom_text(aes(label = lab), vjust = 0, nudge_y = 0.25, size = 2) +
  # geom_line(aes(x = Position, y = 0, group = lineGroup), data = msaDFsnps3) +
  geom_segment(aes(x = PosM1, xend = Position, y = 0, yend = 0), data = msaDFsnps3) +
  xlim(c(1, 3948)) +
  ylim(c(0, 1.5)) +
  facet_wrap(~Sample, ncol = 1, strip.position = "right") +
  theme_void()
ggsave("figures/lollipop_snps.pdf", height = 7, width = 11.5)



ordd <- readr::read_lines("data/lollipop_order.txt")
ordd[ordd == "REF"] <- paste0(ordd[ordd == "REF"], 1:sum(ordd == "REF"))

msaDFsnps4 <- msaDFsnps3[msaDFsnps3$Sample %in% ordd, ]
mmm <- data.frame(Sample = "REF", Position = 1:3947,
    IsAncient = NA, Ref = NA, Alt = NA, y = 0, lab = NA,
    lineGroup = 1, PosM1 = 0:3946)
for (i in 1:sum(grepl("REF", ordd))) {
  mmm_i <- mmm
  mmm_i$Sample <- paste0("REF", i)
  msaDFsnps4 <- rbind(msaDFsnps4, mmm_i)
}
for (i in c("SAMEA6661724", "SAMEA6661722", "SAMEA6661726", "SAMEA5847473")) {
  mmm_i <- mmm
  mmm_i$Sample <- i
  msaDFsnps4 <- rbind(msaDFsnps4, mmm_i)
}
msaDFsnps4$Sample <- factor(msaDFsnps4$Sample, levels = ordd)

ggplot(msaDFsnps4[msaDFsnps4$y > 0, ],
  aes(x = Position, y = y)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = y), lwd = 0.5) +
  geom_point(shape = 21, fill = "white", colour = "black", size = 2) +
  geom_text(aes(label = lab), vjust = 0, nudge_y = 0.25, size = 2) +
  geom_segment(aes(x = PosM1, xend = Position, y = 0, yend = 0), data = msaDFsnps4) +
  xlim(c(1, 3948)) +
  ylim(c(0, 1.5)) +
  facet_wrap(~Sample, ncol = 1, strip.position = "right") +
  theme_void()
ggsave("figures/lollipop_snps2.pdf", height = 13, width = 11.5)

ggplot(msaDFsnps4[msaDFsnps4$y > 0, ],
  aes(x = Position, y = y)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = y), lwd = 0.5) +
  geom_point(shape = 21, fill = "white", colour = "black", size = 2) +
  geom_segment(aes(x = PosM1, xend = Position, y = 0, yend = 0), data = msaDFsnps4) +
  xlim(c(1, 3948)) +
  ylim(c(0, 1.5)) +
  facet_wrap(~Sample, ncol = 1, strip.position = "right") +
  theme_void()
ggsave("figures/lollipop_snps2_nolab.pdf", height = 13, width = 11.5)

###############################################################################


msaAA <- t(apply(msaOriginal, 1,
    function(x) tapply(x, rep(1:(length(x) / 3), each = 3), 
      function(y) paste0(y, collapse = ""))))
replaceNA <- function(x) { x[is.na(x)] <- "-" ; x }
msaAA <- t(apply(msaAA, 1, function(x) replaceNA(Biostrings::GENETIC_CODE[x])))
rownames(msaAA) <- gsub(".consensus.masked", "", rownames(msaAA), fixed = TRUE)

msaAAsub <- msaAA
msaAAref <- msaAA["Reference", ]
for (i in 1:ncol(msaAAsub)) {
  msaAAsub[, i][msaAA[, i] == msaAAref[i]] <- "ref"
  msaAAsub[, i][msaAA[, i] != msaAAref[i] & msaAA[, i] != "-"] <- "sub"
  msaAAsub[, i][msaAA[, i] == "-"] <- "del"
}

rn <- rownames(msaAA)
rnAnc <- rn[grepl("SAM", rn)]
rn <- rn[!rn %in% rnAnc]

for (i in 1:ncol(msaAAsub)) {
  colIsAnc <- msaAAsub[rnAnc, i]
  colIsAnc <- if ("sub" %in% colIsAnc && !"sub" %in% msaAAsub[rn, i]) TRUE else FALSE
  if (colIsAnc) {
    msaAAsub[rnAnc, i][msaAAsub[rnAnc, i] == "sub"] <- "anc"
  }
}
msaAAsubDF <- reshape2::melt(msaAAsub)
msaAAsubDF$Var1 <- as.character(msaAAsubDF$Var1)
msaAAsubDF$Var2 <- as.character(msaAAsubDF$Var2)
msaAAsubDF$residue <- reshape2::melt(msaAA)$value
msaAAsubDF$ref <- rep(msaAAref, each = nrow(msaAAsub))

msaAAsubDFanc <- msaAAsubDF[msaAAsubDF$value == "anc", ]
msaAAsubDFanc$substitution <- paste0(msaAAsubDFanc$ref, msaAAsubDFanc$Var2, msaAAsubDFanc$residue)

table(duplicated(msaAAsubDFanc$substitution[msaAAsubDFanc$Var1 %in% c("SAMEA3486783", "SAMN02727818")]))

# Chinchorro: SAMEA3486783 
# El-Yaral: SAMN02727821
# Chiribaya: SAMN02727818

mmm <- msaDFsnps4[msaDFsnps4$Sample %in% c("SAMEA3486783", "SAMN02727821", "SAMN02727818"), ]
mmm <- mmm[!is.na(mmm$IsAncient), ]

mmm2 <- mmm[mmm$IsAncient == "Yes", ]
mmm3 <- mmm[mmm$IsAncient == "No", ]

nuc2aaPos <- rep(1:ncol(msaAA), each = 3)

mmm2$AApos <- nuc2aaPos[mmm2$Position]

readr::write_tsv(mmm2[order(mmm2$Sample), ], "results/snps_tab.txt")

readr::write_tsv(
  msaAAsubDF[msaAAsubDF$Var1 %in% c("SAMEA3486783", "SAMN02727821", "SAMN02727818"), ] %>%
    .[paste0(.$Var1, .$Var2) %in% paste0(mmm2$Sample, mmm2$AApos), ] %>%
    .[order(.$Var1), ],
  "results/snps_aa_tab.txt")

mmm7 <- mmm[mmm$Sample %in% c("SAMEA3486783", "SAMN02727821", "SAMN02727818"), ]
mmm7$IsShared <- "No"
mmm7$lab <- paste0(mmm7$Ref, mmm7$Position, mmm7$Alt)
mmm7$IsShared[mmm7$lab %in% mmm7$lab[duplicated(mmm7$lab)]] <- "Yes"

mmm7 <- mmm7[order(mmm7$Sample, mmm7$Position), ]
mmm7$AApos <- nuc2aaPos[mmm7$Position]

mmm7 <- mmm7[, c("Sample", "Position", "Ref", "Alt", "AApos", "IsAncient", "IsShared")]

mmm7$RefAA <- msaAA["Reference", mmm7$AApos]
mmm7$AltAA <- mmm7$RefAA

mmm7$AltAA[mmm7$Sample == "SAMEA3486783"] <- msaAA["SAMEA3486783", mmm7$AApos[mmm7$Sample == "SAMEA3486783"]]
mmm7$AltAA[mmm7$Sample == "SAMN02727821"] <- msaAA["SAMN02727821", mmm7$AApos[mmm7$Sample == "SAMN02727821"]]
mmm7$AltAA[mmm7$Sample == "SAMN02727818"] <- msaAA["SAMN02727818", mmm7$AApos[mmm7$Sample == "SAMN02727818"]]

readr::write_tsv(mmm7, "results/snps_select_aa_tab_all.txt")
readr::write_tsv(
  mmm[mmm$Sample %in% c("SAMEA3486783", "SAMN02727821", "SAMN02727818"), ],
  "results/snps_aa_tab_all.txt")


table(duplicated(mmm2$lab[mmm2$Sample %in% c("SAMEA3486783", "SAMN02727821")]))
table(duplicated(mmm2$Position[mmm2$Sample %in% c("SAMEA3486783", "SAMN02727821")]))

table(duplicated(mmm2$lab[mmm2$Sample %in% c("SAMEA3486783", "SAMN02727818")]))
table(duplicated(mmm2$Position[mmm2$Sample %in% c("SAMEA3486783", "SAMN02727818")]))

table(duplicated(mmm2$lab[mmm2$Sample %in% c("SAMN02727821", "SAMN02727818")]))
table(duplicated(mmm2$Position[mmm2$Sample %in% c("SAMN02727821", "SAMN02727818")]))


table(duplicated(mmm3$lab[mmm3$Sample %in% c("SAMEA3486783", "SAMN02727821")]))
table(duplicated(mmm3$Position[mmm3$Sample %in% c("SAMEA3486783", "SAMN02727821")]))

table(duplicated(mmm3$lab[mmm3$Sample %in% c("SAMEA3486783", "SAMN02727818")]))
table(duplicated(mmm3$Position[mmm3$Sample %in% c("SAMEA3486783", "SAMN02727818")]))

table(duplicated(mmm3$lab[mmm3$Sample %in% c("SAMN02727821", "SAMN02727818")]))
table(duplicated(mmm3$Position[mmm3$Sample %in% c("SAMN02727821", "SAMN02727818")]))

#################################################################################

msaAnc <- msa[grepl("^SAM", rownames(msa)), ]

# % coverage of ancient toxin genes
100 - range((apply(msaAnc, 1, \(x) sum(is.na(x))) / ncol(msaAnc)) * 100)

# Number w/o ancient SNPs
NoAnc <- apply(msaAnc, 1, \(x) sum(x[!is.na(x)] == "Yellow"))
sum(NoAnc == 0)

# %ID of those w/ ancient SNPs
range(apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(x[!is.na(x)] == "Black")) / 
  apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(!is.na(x)))) * 100

# # of unique ancient SNPs
sum(apply(msaAnc, 2, \(x) any(x[!is.na(x)] == "Yellow")))

#########################

msaNoAnc <- msa[!grepl("^SAM", rownames(msa)), -ncol(msa)]
100 - range((apply(msaNoAnc, 1, \(x) sum(is.na(x))) / ncol(msaNoAnc)) * 100)
sort(apply(msaNoAnc, 1, \(x) sum(x[!is.na(x)] == "Black")) / 
  apply(msaNoAnc, 1, \(x) sum(!is.na(x)))) * 100

msaAnc <- msa[grepl("^SAM", rownames(msa)), -ncol(msa)]

# % coverage of ancient toxin genes
100 - range((apply(msaAnc, 1, \(x) sum(is.na(x))) / ncol(msaAnc)) * 100)
a1 <- 100 - unname((apply(msaAnc, 1, \(x) sum(is.na(x))) / ncol(msaAnc)) * 100)

# Number w/o ancient SNPs
NoAnc <- apply(msaAnc, 1, \(x) sum(x[!is.na(x)] == "Yellow"))
sum(NoAnc == 0)

# %ID of those w/ ancient SNPs
range(apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(x[!is.na(x)] == "Black")) / 
  apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(!is.na(x)))) * 100

a1 <- 100 - unname((apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(is.na(x))) / ncol(msaAnc)) * 100)
a2 <- unname(apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(x[!is.na(x)] == "Black")) / 
  apply(msaAnc[NoAnc > 0, ], 1, \(x) sum(!is.na(x)))) * 100
a3 <- data.frame(Pct.Cov = a1, Pct.ID = a2)
a3 <- a3[order(a3$Pct.ID), ]

b1 <- 100 - unname((apply(msaAnc[NoAnc == 0, ], 1, \(x) sum(is.na(x))) / ncol(msaAnc)) * 100)
b2 <- unname(apply(msaAnc[NoAnc == 0, ], 1, \(x) sum(x[!is.na(x)] == "Black")) / 
  apply(msaAnc[NoAnc == 0, ], 1, \(x) sum(!is.na(x)))) * 100
b3 <- data.frame(Pct.Cov = b1, Pct.ID = b2)
b3 <- b3[order(b3$Pct.ID), ]

# # of unique ancient SNPs
sum(apply(msaAnc, 2, \(x) any(x[!is.na(x)] == "Yellow")))

msa9 <- msaOriginal[rownames(msaAnc), -ncol(msaOriginal)]
msa8 <- msa9
for (i in 1:ncol(msa8)) {
  for (j in 1:nrow(msa8)) {
    if (is.na(msaAnc[j, i])) {
      msa8[j, i] <- "."
    }  else if (msaAnc[j, i] != "Yellow") {
      msa8[j, i] <- "-"
    }
  }
}

msa8 <- msa8[!apply(msa8, 1, \(x) all(x %in% c("-", "."))), ]
msa8 <- msa8[, !apply(msa8, 2, \(x) all(x %in% c("-", ".")))]

msa7 <- apply(msa8, 1, \(x) paste0(x, collapse = ""))
msa7 <- msa7[!duplicated(msa7)]

Atab <- readr::read_tsv("data/TeNT_pctCov.txt")

rownames(msaAnc) <- gsub(".consensus.masked", "", rownames(msaAnc))
for (i in 1:20) {
  tmp <- msaAnc[Atab$BioSample[i], ]
  Atab[[4]][i] <- sum(tmp[!is.na(tmp)] == "Black") / sum(!is.na(tmp))
}

Atab[[2]] <- round(Atab[[2]] * 100, 2)
Atab[[4]] <- round(Atab[[4]] * 100, 2)
readr::write_tsv(Atab, "data/TeNT_pctID.txt")

msaO <- msaOriginal[, -ncol(msaOriginal)]
rownames(msaO) <- gsub(".consensus.masked", "", rownames(msaO))
msaO1 <- msaO[grepl("SAM", rownames(msaO)), ]
msaO2 <- msaO[!grepl("SAM", rownames(msaO)), ]

ref <- msaO["Reference", ]
refID <- numeric(nrow(msaO))
for (i in 1:nrow(msaO)) {
  tmp <- msaO[i, ]
  gaps <- is.na(tmp) | tmp == "-"
  tmp_ref <- ref[!gaps]
  tmp <- tmp[!gaps]
  refID[i] <- sum(tmp_ref == tmp) / length(tmp)
}
names(refID) <- rownames(msaO)

Atab[[4]] <- round(refID[Atab$BioSample] * 100, 2)
readr::write_tsv(Atab, "data/TeNT_pctID.txt")

mmm <- msaO[c("SAMEA3486783", "SAMN02727821"), ]
mmm <- mmm[, apply(mmm, 2, \(x) !any(is.na(x) | x == "-"))]
sum(mmm[1, ] == mmm[2, ]) / ncol(mmm)

cmp_f <- function(x, y) {
  gaps <- apply(matrix(c(x, y), ncol = 2), 1, \(z) any(is.na(z) | z == "-"))
  if (all(gaps)) {
    message("Nothing to compare")
    0
  } else {
    sum(x[!gaps] == y[!gaps]) / sum(!gaps)
  }
}

toCmp <- readr::read_lines("data/TeNT_75PctID.txt")
msaO1 <- msaO[toCmp, ]

msaO_cmp <- matrix(0, ncol = nrow(msaO2), nrow = nrow(msaO1),
  dimnames = list(rownames(msaO1), rownames(msaO2)))
for (i in 1:ncol(msaO_cmp)) {
  for (j in 1:nrow(msaO_cmp)) {
    msaO_cmp[j, i] <- cmp_f(
      msaO[colnames(msaO_cmp)[i], ],
      msaO[rownames(msaO_cmp)[j], ])
  }
}

msaO_cmp_W <- data.frame(
  Sequence = rownames(msaO_cmp),
  HighestPctID = NA_real_,
  HighestPctIDSequence = NA_character_
)
for (i in 1:nrow(msaO_cmp_W)) {
  msaO_cmp_W$HighestPctID[i] <- msaO_cmp[i, which.max(msaO_cmp[i, ])]
  msaO_cmp_W$HighestPctIDSequence[i] <- paste0(
    colnames(msaO_cmp)[msaO_cmp[i, ] == msaO_cmp_W$HighestPctID[i]],
    collapse = ", "
  )
}
msaO_cmp_W$HighestPctID <- round(msaO_cmp_W$HighestPctID * 100, 2)
readr::write_tsv(msaO_cmp_W, "results/TeNT_75PctID_NearestMatch.tsv")

mmz <- msaO[c("SAMN02799091", "SAMN02799089"), ]
mmz <- apply(mmz, 1, \(x) paste0(x, collapse = ""))
Biostrings::writeXStringSet(Biostrings::DNAStringSet(mmz), "results/LowPctIDSeqs.fa")

