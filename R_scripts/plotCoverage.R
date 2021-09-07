#!/usr/bin/env Rscript
# Benjamin Jean-Marie Tremblay
# June 22nd, 2021
# version 8
# 
#
# ChangeLog:
#
# June 22nd, 2021 (version 8)
# - Fix weird --names position bug (with terrible hack!)
#
# June 21st, 2021 (version 7)
# - Try and fix ylim values of --names
# - Add a filename check to --track-order
#
# June 20th, 2021 (version 6)
# - New --names option: plot track names
# - New --names-cex option: track name size
# - New --open-angle: opening angle for plotting track names
# - New --track-order: change ordering of tracks
#
# May 26th, 2021 (version 5)
# - Fixed --colour arg parsing picking up --hl-colours arg
#
# May 26th, 2021 (version 4)
# - New --hl-features option: colour specific genes/features
# - New --hl-colours: customize the --hl-features colours
#
# May 14th, 2021 (version 3)
# - New --R-lib option: where to install/load packages from
# - New --plot-sum option: whether to plot an outer ring of total summed reads
# - And also --sum-colour, --sum-height
# - New --no-border: hide area plot borders
#
# April 15th, 2021 (version 2)
# - Typo in usage
# - Check for empty BAM list
# - Tidy help check regex
# - Only read .bam files in --bam folder

message("Running plotCoverage.R v8")
message("-------------------------")
message("")

print_usage <- function() {
  message("Usage:    plotCoverage --gff=file.gff --bam=bamFolder/ [optional args]")
  message("")
  message("  --gff           Name of GFF file. (required)")
  message("  --bam           Name of folder containing .bam files. (required)")
  message("  --out           Output filename. (default: --out=plot.pdf)")
  message("  --plot-genes    Whether to plot genes. (default: --plot-genes=true)")
  message("  --gene-fname    Feature name to plot if --plot-genes=true.")
  message("                  (default: --gene-fname=gene)")
  message("  --names         Plot the names of sequences/tracks by opening up the circos")
  message("                  and plotting the names besides the ends of the sequences.")
  message("                  If --track-order is not provided, then the file names are")
  message("                  used. Note this is only intended for single contig or")
  message("                  chromosome-containing bam files. (default: --names=false)")
  message("  --names-cex     Sequence/track name cex. (default: --names-cex=0.25)")
  message("  --open-angle    The angle of the opening when --names=true. (default:")
  message("                  --open-angle=45)")
  message("  --track-order   Path to a tsv file, with the first column being the")
  message("                  bam file name and the second column being the track name.")
  message("                  The order of rows will be used as the order of tracks.")
  message("                  The tsv file must not have column names.")
  message("                  sequence/track. Used to control ordering of tracks.")
  message("  --hl-features   Highlight specific genes if --plot-genes=true. Provide")
  message("                  comma-separated coordinates in the form sequence:start..stop")
  message("                  (e.g. --hl-features=Chr1:147..752,Chr2:6260..7162)")
  message("  --hl-colours    Comma-separated list of colours to be used alongside")
  message("                  --hl-features. If only a single colour is provided, then")
  message("                  it will be recycled for all highlighted features; otherwise")
  message("                  the number of colours must match the number of features.")
  message("                  (default: --hl-colours=red)")
  message("  --gff-source    GFF data source to be used for getting chromosome/contig")
  message("                  size and gene coordinates. (default: --gff-source=RefSeq)")
  message("  --bin-size      Bin size for calculating mean read coverage.")
  message("                  (default: --bin-size=1000)")
  message("  --max-pctile    Max read count percentile to plot. Helps with unbalanced")
  message("                  coverage. (default: --max-pctile=90)")
  message("  --cov-height    Height of coverage plots, as a fraction of the total radius.")
  message("                  (default: 1 divided by double the number of BAM files)")
  message("  --gene-height   Height of gene plot, as a fraction of the total radius.")
  message("                  (default: --gene-height=0.01)")
  message("  --single-ylim   Whether to have all coverage plots share the same ylim.")
  message("                  (default: --single-ylim=true)")
  message("  --label-cex     Chromosome/contig text cex. (default: --label-cex=0.4)")
  message("  --tick-cex      Tick label text cex. (default: --label-cex=0.3)")
  message("  --abbr-labels   Abbreviate chromosome/contig names.")
  message("                  (default: --abbr-labels=true)")
  message("  --min-size      Minimum chromosome/contig size. (default: --min-size=10000)")
  message("  --tick-size     Tick interval size. (default: --tick-size=50000)")
  message("  --colour        Area plot fill colour. (default: --colour=blue)")
  message("  --no-border     Don't plot area plot borders. (default: --no-border=false)")
  message("  --plot-sum      Plot the sum of reads across all samples as an outer ring.")
  message("                  (default: --plot-sum=false)")
  message("  --sum-height    Height of --plot-sum track, as a fraction the total radius.")
  message("                  (default: --sum-height=0.04)")
  message("  --sum-colour    Area plot fill colour for --plot-sum outer ring. (default:")
  message("                  --sum-colour=darkgreen)")
  message("  --R-lib         Location of the R library to install and load packages")
  message("                  to/from. (Not setting this will use the default location.)")
  message("")
}

args <- commandArgs(TRUE)
check_for_help <- function(x) {
  grepl("^--help", x) || grepl("^-help", x) || grepl("^--h", x) || grepl("^-h", x)
}
if (!length(args) || check_for_help(args)) {
  if (!interactive()) {
    print_usage()
    q("no")
  }
}

if (!any(grepl("--gff=", args)))
  stop("Missing GFF file.")
if (!any(grepl("--bam=", args)))
  stop("Missing BAM folder.")

get_arg <- function(x) {
  gsub(paste0("^--", x, "="), "", args[grepl(paste0(x, "="), args)])
}

get_arg_opt <- function(x, def) {
  if (any(grepl(paste0("^--", x, "="), args)))
    get_arg(x)
  else
    def
}

argList <- list(
  gff = get_arg("gff"),
  bam = get_arg("bam"),
  out = get_arg_opt("out", "plot.pdf"),
  plot_genes = as.logical(toupper(get_arg_opt("plot-genes", "true"))),
  hl_features = suppressMessages(get_arg_opt("hl-features", NULL)),
  hl_colours = suppressMessages(get_arg_opt("hl-colours", "red")),
  gff_source = get_arg_opt("gff-source", "RefSeq"),
  bin_size = as.integer(get_arg_opt("bin-size", "1000")),
  max_pctile = as.numeric(get_arg_opt("max-pctile", "90")) / 100,
  label_cex = as.numeric(get_arg_opt("label-cex", "0.4")),
  tick_cex = as.numeric(get_arg_opt("tick-cex", "0.3")),
  abbr_labels = as.logical(toupper(get_arg_opt("abbr-labels", "true"))),
  min_size = as.integer(get_arg_opt("min-size", "10000")),
  tick_size = as.integer(get_arg_opt("tick-size", "50000")),
  colour = get_arg_opt("colour", "blue"),
  gene_fname = get_arg_opt("gene-fname", "gene"),
  cov_height = suppressMessages(get_arg_opt("cov-height", NULL)),
  single_ylim = as.logical(toupper(get_arg_opt("single-ylim", "true"))),
  gene_height = as.numeric(get_arg_opt("gene-height", "0.01")),
  plot_sum = as.logical(toupper(get_arg_opt("plot-sum", "false"))),
  names = as.logical(toupper(get_arg_opt("names", "false"))),
  names_cex = as.numeric(get_arg_opt("names-cex", "0.25")),
  open_angle = as.numeric(get_arg_opt("open-angle", "45")),
  track_order = suppressMessages(get_arg_opt("track-order", NULL)),
  sum_height = as.numeric(get_arg_opt("sum-height", "0.04")),
  sum_colour = get_arg_opt("sum-colour", "darkgreen"),
  R_lib = suppressMessages(get_arg_opt("R-lib", NULL)),
  no_border = as.logical(toupper(get_arg_opt("no-border", "false")))
)

if (!nchar(argList$gff)) stop("Incorrect --gff")
if (!nchar(argList$bam)) stop("Incorrect --bam")
if (!nchar(argList$out)) stop("Incorrect --out")
if (is.na(argList$plot_genes)) stop("Incorrect --plot-genes")
if (is.na(argList$plot_sum)) stop("Incorrect --plot-sum")
if (!nchar(argList$gff_source)) stop("Incorrect --gff-source")
if (is.na(argList$bin_size)) stop("Incorrect --bin-size")
if (is.na(argList$max_pctile)) stop("Incorrect --max-pctile")
if (is.na(argList$label_cex)) stop("Incorrect --label-cex")
if (is.na(argList$tick_cex)) stop("Incorrect --tick-cex")
if (is.na(argList$abbr_labels)) stop("Incorrect --abbr-labels")
if (is.na(argList$min_size)) stop("Incorrect --min-size")
if (is.na(argList$tick_size)) stop("Incorrect --tick-size")
if (!nchar(argList$colour)) stop("Incorrect --colour")
if (!nchar(argList$hl_colours)) stop("Incorrect --hl-colours")
if (!nchar(argList$sum_colour)) stop("Incorrect --sum-colour")
if (is.na(argList$single_ylim)) stop("Incorrect --single_ylim")
if (is.na(argList$gene_height)) stop("Incorrect --gene-height")
if (is.na(argList$sum_height)) stop("Incorrect --sum-height")
if (is.na(argList$no_border)) stop("Incorrect --no-border")
if (is.na(argList$names)) stop("Incorrect --names")
if (is.na(argList$names_cex)) stop("Incorrect --names-cex")
if (is.na(argList$open_angle)) stop("Incorrect --open-angle")

if (!is.null(argList$track_order) &&
    (is.na(argList$track_order) || !nchar(argList$track_order))) {
  stop("Incorrect --track-order")
}

if (!is.null(argList$R_lib) && is.na(argList$R_lib)) {
  stop("Incorrect --R-lib")
} else if (!is.null(argList$R_lib)) {
  lpaths <- .libPaths()
  .libPaths(c(argList$R_lib, lpaths))
}

if (!is.null(argList$cov_height) && is.na(argList$cov_height)) {
  stop("Incorrect --cov-height")
} else if (!is.null(argList$cov_height)) {
  argList$cov_height <- as.numeric(argList$cov_height)
}

if (!is.null(argList$hl_features) && is.na(argList$hl_features)) {
  stop("Incorrect --hl-features")
}

message("Loading packages...")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager.")
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

check_pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    message("Installing missing dependencies.")
    BiocManager::install(x)
  }
}

check_pkg("circlize")
check_pkg("Rsamtools")
check_pkg("GenomicFeatures")
check_pkg("readr")

suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicFeatures))

message("Loading GFF...")

gff <- suppressMessages(readr::read_tsv(argList$gff, col_names = FALSE, comment = "#"))

gff2 <- gff[gff[[2]] == argList$gff_source, ]
if (!nrow(gff2))
  stop("Couldn't find any rows in GFF with source \"", argList$gff_source, "\"")
gff2 <- gff2[, c(1, 4:5, 3)]
colnames(gff2) <- c("chr", "start", "end", "feature")
gff2 <- gff2[!is.na(gff2$chr), ]

contigs <- gff2[gff2$feature == "region", ]
if (!nrow(contigs))
  stop("Couldn't find any rows in GFF with feature \"region\"")
contigsGR <- GRanges(contigs)[, -1]
seqlengths(contigsGR) <- end(contigsGR)

genes <- gff2[gff2$feature == argList$gene_fname, ]
if (!nrow(genes) && argList$plot_genes)
  stop("Couldn't find any rows in GFF with feature \"", argList$gene_fname, "\"")

# This is the bin width, change as desided
contigBins <- tileGenome(seqlengths(contigsGR), tilewidth = argList$bin_size,
  cut.last.tile.in.chrom = TRUE)

# Get BAM coverage from files
#   meanMaxQuant: controls the max percentile for coverage)
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
  out$MeanReads[out$MeanReads > quant] <- quant

  out <- rbind(out[1, ], out)
  out$start[1] <- 0
  out$end[1] <- 1
  out$MeanReads[1] <- 0

  out

}

message("Loading BAM files...")

# Load bam files into a list (just one for now)
bamfiles <- list.files(argList$bam, "*.bam$", full.names = TRUE)
bamfilesShort <- list.files(argList$bam, "*.bam$", full.names = FALSE)
if (!length(bamfiles))
  stop("Couldn't find any BAM files in folder \"", argList$bam, "\"")
bam <- lapply(bamfiles, function(x) get_bam_dat(x, meanMaxQuant = argList$max_pctile))

if (argList$single_ylim) {
  max_ylims <- vapply(bam, function(x) max(x[[4]]), numeric(1))
  for (i in seq_along(bam)) {
    bam[[i]][[4]] <- max(max_ylims) - bam[[i]][[4]]
    bam[[i]][[4]] <- bam[[i]][[4]] / max(max_ylims)
  }
} else {
  for (i in seq_along(bam)) {
    bam[[i]][[4]] <- max(bam[[i]][[4]]) - bam[[i]][[4]]
    bam[[i]][[4]] <- bam[[i]][[4]] / max(bam[[i]][[4]])
  }
}

# Nicer labels
if (argList$abbr_labels) {
  contig2num <- structure(seq_len(nrow(contigs)), names = contigs$chr)
} else {
  contig2num <- structure(contigs$chr, names = contigs$chr)
}

# Track height for coverage plots
if (is.null(argList$cov_height)) {
  cov_height <- 1 / (length(bam) * 2)
} else {
  cov_height <- argList$cov_height
}

# Sum of reads across all samples
if (argList$plot_sum) {
  bam2 <- as.data.frame(lapply(bam, function(x) x[["MeanReads"]]))
  bam2 <- rowSums(bam2)
  bamSum <- bam[[1]]
  bamSum$MeanReads <- bam2
  bamSum$MeanReads <- max(bamSum$MeanReads) - bamSum$MeanReads
}

genes$colour <- "black"

if (argList$plot_genes && !is.null(argList$hl_features)) {
  hl_f <- strsplit(argList$hl_features, ",", TRUE)[[1]]
  hl_c <- strsplit(argList$hl_colours, ",", TRUE)[[1]]
  hl_c <- rep_len(hl_c, length(hl_f))
  new_g <- data.frame(
    chr = rep(NA_character_, length(hl_f)),
    start = rep(NA_real_, length(hl_f)),
    end = rep(NA_real_, length(hl_f)),
    feature = argList$gene_fname,
    colour = unlist(hl_c)
  )
  for (i in seq_along(hl_f)) {
    hl_f_tmp <- strsplit(hl_f[i], ":", TRUE)[[1]]
    hl_f_chr <- hl_f_tmp[1]
    hl_f_crd <- as.numeric(strsplit(hl_f_tmp[2], "..", TRUE)[[1]])
    new_g$chr[i] <- hl_f_chr
    new_g$start[i] <- hl_f_crd[1]
    new_g$end[i] <- hl_f_crd[2]
  }
  genes <- rbind(genes, new_g)
}

if (argList$names) {
  if (!is.null(argList$track_order)) {
    names(bam) <- bamfilesShort
    tnames <- suppressMessages(readr::read_tsv(argList$track_order, col_names = FALSE))
    err <- FALSE
    if (any(!basename(tnames[[1]]) %in% names(bam))) {
      message("Error: couldn't find the following BAM files from the --track-order file:")
      colnames(tnames) <- c("FileName", "NewName")
      print(as.data.frame(tnames[!basename(tnames[[1]]) %in% names(bam), ]))
      err <- TRUE
    }
    if (any(!names(bam) %in% basename(tnames[[1]]))) {
      message("Error: the following BAM filenames couldn't be found in the --track-order file:")
      print(names(bam)[!names(bam) %in% basename(tnames[[1]])])
      err <- TRUE
    }
    if (err) q()
    bam <- bam[basename(tnames[[1]])]
    names(bam) <- tnames[[2]]
  } else {
    names(bam) <- gsub(".bam$", "", bamfilesShort)
  }
}

# circlize --------------------------------------------------------------------

message("Making plot...")

circos.clear()
pdf(argList$out)

circos.par(
  start.degree = if (argList$names) 90 - argList$open_angle else 90,
  cell.padding = c(0, 0, 0, 0),
  gap.degree = if (argList$names) argList$open_angle else 0.8,
  points.overflow.warning = FALSE,
  track.margin = c(0.0075, 0)
)

ideogram <- cbind(contigs[, 1:3], "A", "gneg")

# Outer axis
circos.initializeWithIdeogram(ideogram,
  plotType = "axis",
  chromosome.index = contigs$chr[contigs$end >= argList$min_size],
  tickLabelsStartFromZero = FALSE, major.by = argList$tick_size,
  axis.labels.cex = argList$tick_cex,
  labels.cex = argList$label_cex)

# Labels
circos.genomicTrack(ideogram,
  ylim = c(0, 1), track.height = 0.001, bg.col = NA, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.text(mean(CELL_META$xlim), 70, contig2num[CELL_META$sector.index],
      niceFacing = TRUE, cex = argList$label_cex)
})

if (argList$plot_sum) {
  # Sum of all BAMs
  circos.genomicTrack(bamSum, track.height = argList$sum_height,
    bg.lwd = .4, bg.border = if (argList$no_border) NA else "black",
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, area = TRUE, baseline = "bottom",
        col = argList$sum_colour, lwd = 0.5, border = NA)
  })
}

if (argList$plot_genes) {
# Genes
  circos.genomicTrackPlotRegion(genes,
    ylim = c(0, 1), bg.border = NA, track.height = argList$gene_height,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, ytop = 1, ytbottom = 0, col = value$colour, border = NA)
  })
}

# Iterate through all the bam files and plot coverage
for (i in seq_along(bam)) {
  circos.genomicTrack(bam[[i]], track.height = cov_height,
    bg.lwd = .4, bg.border = if (argList$no_border) NA else "black",
    panel.fun = function(region, value, ...) {
      if (argList$names) {
        x_tmp <- ideogram[[3]][ideogram[[1]] == CELL_META$sector.index]
        # y=.5 SHOULD be the proper setting, but for Harold's BAM files it
        # bugs out completely. I have zero idea why. Setting y=1 and
        # adj=c(0,1.1) is a hacky way of getting around it.
        circos.text(x_tmp * 1.00135, y = 1, #y = .5,
          labels = names(bam)[i],
          # cex = argList$names_cex, facing = "downward", adj = c(0, 0.5))
          cex = argList$names_cex, facing = "downward", adj = c(0, 1.1))
      }
      circos.genomicLines(region, value, area = TRUE, baseline = "top",
        col = argList$colour[1], lwd = 0.5, border = NA)
  })
}

dev.off()

message("Done.")
