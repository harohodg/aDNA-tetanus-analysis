# Benjamin Jean-Marie Tremblay
# 2021-10-31

library(ggplot2)

msa <- Biostrings::readDNAStringSet("data/TeNT-MSA-nucl.mfa.longer.fasta")
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

rowOrder <- c(
  "SAMEA6661724",
  "2017.061",
  "SAMEA6661722",
  "SAMEA6661726",
  "512-15",
  "SAMEA5847473",
  "202.15",
  "SAMEA3937653",
  "SAMD00041000",
  "SAMD00041001",
  "SAMN05991104",
  "SAMEA104281224",
  "778.17",
  "SAMEA5847472",
  "SAMN02727821",
  "SAMEA3486783",
  "SAMN02727818",
  "ATCC_453",
  "1337",
  "SAMEA6490841",
  "63.05",
  "TMB2",
  "SAMEA5847426",
  "132CV",
  "12124569",
  "SAMEA104281220",
  "C2",
  "SAMEA104281226",
  "SAMEA5847432",
  "ATCC_9441",
  "SAMEA104281221",
  "SAMEA104441581"
)

rowOrder <- rownames(msa)[pmatch(rowOrder, rownames(msa))]
msa <- msa[rowOrder, ]

SNPcount <- apply(msa, 1, \(x) sum(x == "Yellow", na.rm = TRUE))
SNPcount <- data.frame(
  BioSample = names(SNPcount),
  AncientSNPsCount = unname(SNPcount)
)
# readr::write_tsv(SNPcount, "ToxinAncientSNPsCount.tsv")

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

# readr::write_tsv(msaDFsnps, "ToxinSNPs.tsv")

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
ggsave("ToxinSNPs.png",
  plot = p, width = 6.75, height = 4.2)
ggsave("ToxinSNPs.pdf",
  plot = p, width = 6.75, height = 4.2)

