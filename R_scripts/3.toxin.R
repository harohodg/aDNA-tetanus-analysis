# Benjamin Jean-Marie Tremblay
# 2021-07-10

library(ggplot2)

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

p <- ggplot(msaDF, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0),
    breaks = c(1000, 2000, 3000),
    labels = c("1KB", "2KB", "3KB")) +
  scale_y_discrete(position = "right") +
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
ggsave("figures/ToxinSNPs.png",
  plot = p, width = 6.75, height = 4.2)
ggsave("figures/ToxinSNPs.pdf",
  plot = p, width = 6.75, height = 4.2)

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
