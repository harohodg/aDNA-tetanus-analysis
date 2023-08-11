library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)
library(ape)

gffF <- "data/E88_plasmid.gff3"
bams <- list.files("data/bam-plasmid/", "*.bam$", full = TRUE)
names(bams) <- gsub("_.*", "", basename(bams))
toxin <- GRanges("NC_004565.1", IRanges(68640, 72587))

get_bam_dat <- function(x, meanMaxQuant = .9) {

  bam <- as.data.frame(scanBam(x, param = ScanBamParam(what = c("rname", "pos", "qwidth"))))
  colnames(bam) <- c("seqname", "start", "width")
  # Drop unused levels carried over in factor column
  bam[[1]] <- as.character(bam[[1]])
  bam$end <- bam$start + bam$width - 1
  bam <- GRanges(bam)

  bam <- binnedAverage(contigBins, coverage(bam), "MeanReads")

  out <- data.frame(
    chr = seqnames(bam),
    start = start(bam),
    end = end(bam)
  )

  out$MeanReads <- bam$MeanReads
  quant <- quantile(out$MeanReads, meanMaxQuant)
  if (quant == 0 && any(out$MeanReads > 0)) quant <- min(out$MeanReads[out$MeanReads > 0])
  out$MeanReads[out$MeanReads > quant] <- quant

  out <- rbind(out[1, ], out)
  out$start[1] <- 0
  out$end[1] <- 1
  out$MeanReads[1] <- 0

  out

}

supp <- readr::read_tsv("data/SuppTable2.txt")
tree3 <- read.tree("data/final-snippy-fasttree-oct23.newick")
ord <- data.frame(row.names = NULL, Sample = NA_character_, Name = tree3$tip.label)
ss <- data.frame(row.names = NULL, id = supp[[1]], name = supp[[4]])
ss$name[14] <- "'Yámana-Tooth'"
ss$name[15] <- "'Vác-Mummy-Tissue'"
ord$Sample <- structure(ss$id, names = ss$name)[ord$Name]

ord$Final <- ord$Sample
ord$Final[is.na(ord$Sample)] <- paste0("Ref", 1:sum(is.na(ord$Sample)))

gff0 <- import(gffF)
gff <- suppressMessages(readr::read_tsv(gffF, col_names = FALSE, comment = "#"))
gff2 <- gff[gff[[2]] == "RefSeq", ]
gff2 <- gff2[, c(1, 4:5, 3)]
colnames(gff2) <- c("chr", "start", "end", "feature")
gff2 <- gff2[!is.na(gff2$chr), ]

contigs <- gff2[gff2$feature == "region", ]
contigsGR <- GRanges(contigs)[, -1]
seqlengths(contigsGR) <- end(contigsGR)

add_ref <- function(x) {
  data.frame(row.names = NULL, chr = "NC_004565.1", start = 1, end = 74082,
    MeanReads = 0, sample = x, toxin = FALSE, colT = FALSE, repA = FALSE)
}
add_all_refs <- function(x) {
  do.call(rbind, lapply(x, add_ref))
}

genes <- gff2[gff2$feature == "gene", ]
genesOld <- genes
genes <- as.data.frame(reduce(GRanges(genes)))

contigBins <- tileGenome(seqlengths(contigsGR), tilewidth = 100,
  cut.last.tile.in.chrom = TRUE)

bamDat <- lapply(bams, function(x) get_bam_dat(x, meanMaxQuant = 0.8))

add_genes <- function() {
  x <- data.frame(row.names = NULL, chr = "NC_004565.1",
    start = c(genes$start, genes$start, genes$start, genes$end, genes$end),
    end = c(genes$start, genes$start, genes$end, genes$end, genes$end),
    MeanReads = c(rep(0, nrow(genes)), rep(1, nrow(genes)*3), rep(0, nrow(genes))),
    sample = "Genes", toxin = FALSE, colT = FALSE, repA = FALSE)
  x$toxin[x$start >= 68640 & x$end <= 72587] <- TRUE
  x$colT[x$start >= 39438 & x$end <= 42413] <- TRUE
  x$repA[x$start >= 20250 & x$end <= 21941] <- TRUE
  x <- x[order(x$start, x$end), ]
  x
}

bamDat2 <- bamDat
for (i in seq_along(bamDat2)) bamDat2[[i]]$sample <- names(bamDat2)[i]
bamDat2 <- do.call(rbind, bamDat2)
rownames(bamDat2) <- NULL
bamDat2$toxin <- FALSE
bamDat2$toxin[bamDat2$start >= 68640 & bamDat2$end <= 72587] <- TRUE
bamDat2$colT <- FALSE
bamDat2$colT[bamDat2$start >= 39438 & bamDat2$end <= 42413] <- TRUE
bamDat2$repA <- FALSE
bamDat2$repA[bamDat2$start >= 20250 & bamDat2$end <= 21941] <- TRUE

bamDat3 <- rbind(add_all_refs(ord$Final[grepl("^Ref\\d+", ord$Final)]), bamDat2)
bamDat3 <- bamDat3[bamDat3$sample %in% ord$Final, ]
bamDat3 <- rbind(add_genes(), bamDat3)
bamDat3$sample2 <- factor(structure(c("Genes", ord$Name),
    names = c("Genes", ord$Final))[bamDat3$sample], levels = c("Genes", ord$Name))

p <- ggplot(bamDat3, aes(x = (end + start) / 2, y = MeanReads)) +
  geom_area() +
  geom_area(data = bamDat3[bamDat3$toxin, ], aes(x = (end + start) / 2, y = MeanReads),
    fill = "red", colour = NA) +
  geom_area(data = bamDat3[bamDat3$colT, ], aes(x = (end + start) / 2, y = MeanReads),
    fill = "blue", colour = NA) +
  geom_area(data = bamDat3[bamDat3$repA, ], aes(x = (end + start) / 2, y = MeanReads),
    fill = "purple", colour = NA) +
  facet_wrap(~sample2, ncol = 1, strip.position = "right", scales = "free_y") +
  theme_void() +
  scale_colour_manual(values = c("black", "red")) +
  theme(
    panel.spacing = unit(0.1, "mm"),
    strip.text.y.right = element_text(hjust = 0, size = unit(6, "pt"))
  )
ggsave("figures/TreeOrder_Colour_FreeYaxis_quantile80_window100bp_nov7.pdf", p,
  width = 9.12, height = 6.12)

