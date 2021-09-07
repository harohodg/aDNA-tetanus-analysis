# Benjamin Jean-Marie Tremblay
# 2021-07-10

msa1 <- Biostrings::readDNAStringSet("data/ANDREW2-consensus_variants.E88_aligned.noInsertions.fastANDREW-a")
msa2 <- Biostrings::readDNAStringSet("data/consensus_variants.E88_aligned.noInsertions.fastANDREW-a")

msa1 <- as.character(msa1)
msa2 <- as.character(msa2)

msa1 <- t(sapply(msa1, \(.) strsplit(., "", TRUE)[[1]]))
msa2 <- t(sapply(msa2, \(.) strsplit(., "", TRUE)[[1]]))
msa2 <- msa2[rownames(msa1), ]

n <- 0

for (i in 1:ncol(msa1)) {
  for (j in 1:nrow(msa1)) {
    p1 <- msa1[j, i]
    p2 <- msa2[j, i]
    if (p1 == "-" || p2 == "-") next
    if (p1 != p2) n <- n + 1
  }
}

msa3 <- Biostrings::readDNAStringSet("data/core.full.cleaned.aln")
msa3 <- as.character(msa3)
msa3 <- t(sapply(msa3, \(.) strsplit(., "", TRUE)[[1]]))
msa3 <- msa3[, 2867891:2871838]
msa3 <- msa3[, ncol(msa3):1]

RC <- c("-" = "-", "A" = "T", "T" = "A", "C" = "G", "G" = "C")
msa3 <- apply(msa3, 1:2, \(x) RC[x])

msa1n <- gsub(".consensus.masked", "", rownames(msa1))
msa1n[msa1n == "E88_Reference"] <- "Reference"
rownames(msa1) <- msa1n

msa3 <- msa3[msa1n, ]

n <- 0
Di <- c()
Dj <- c()

for (i in 1:ncol(msa1)) {
  for (j in 1:nrow(msa3)) {
    p1 <- msa1[j, i]
    p2 <- msa3[j, i]
    if (p1 == "-" || p2 == "-") next
    if (p1 != p2) {
      n <- n + 1
      Di <- c(Di, i)
      Dj <- c(Dj, j)
    }
  }
}

Diffs <- data.frame(
  Sample = rownames(msa3)[Dj],
  Pos = Di,
  OldBase = NA_character_,
  NewBase = NA_character_
)
for (i in 1:nrow(Diffs)) {
  Diffs$OldBase[i] <- msa3[Dj[i], Di[i]]
  Diffs$NewBase[i] <- msa1[Dj[i], Di[i]]
}
readr::write_tsv(Diffs, "results/CoreFullCleaned_vs_MikeMSA_Diffs.tsv")
