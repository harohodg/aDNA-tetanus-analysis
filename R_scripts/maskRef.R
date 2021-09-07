# Benjamin Jean-Marie Tremblay
# 2021-05-24

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)

# Read in BAM file
sp <- ScanBamParam(what = c("rname", "strand", "pos", "qwidth"))
bam <- scanBam("./Archive/2017.061.bowtie.sorted.bam", param = sp)[[1]]
bam <- as.data.frame(bam)

colnames(bam) <- c("seqnames", "strand", "start", "width")

# To GRanges
bam$stop <- bam$start + bam$width - 1
bam$strand <- "+"
bam <- GRanges(bam)

# Read in reference
seqs <- readDNAStringSet("./Archive/GCF_004115495.1_ASM411549v1_genomic.fna")

# In this case the names don't quite match
names(seqs) <- vapply(strsplit(names(seqs), " ", TRUE), function(x) x[1], character(1))

# Now format the two objects to have the same metadata
seqlengths(seqs) <- structure(width(seqs), names = names(seqlengths(seqs)))
seqlengths(bam) <- seqlengths(seqs)

# Merge overlapping ranges
bam <- sortSeqlevels(bam)
bam <- sort(bam)
bam <- reduce(bam)

# Some out-of-bound ranges detected. Not sure why as I'm not familiar with this
# data, but I will just trim them off.
bam <- trim(bam)

# What we want to keep
seqs.good <- seqs[bam]

# And generate our bunk sequences
new.seqs <- DNAStringSet(vapply(width(seqs),
    function(x) paste0(rep("N", x), collapse = ""), character(1)))
names(new.seqs) <- names(seqs)
seqlengths(new.seqs) <- seqlengths(seqs)

# Replace our good ranges
bam.split <- split(bam, seqnames(bam))
seqs.good.split <- split(seqs.good, names(seqs.good))
for (i in seq_along(bam.split)) {
  seq.name <- names(bam.split)[i]
  new.seqs[[seq.name]] <- replaceAt(new.seqs[[seq.name]],
    IRanges(start(bam.split[[i]]), width = width(bam.split[[i]])),
    as.character(seqs.good.split[[seq.name]]))
}

# Done!
writeXStringSet(new.seqs, "masked_seqs.fna")
