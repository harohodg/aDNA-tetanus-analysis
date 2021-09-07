# Ben's scripts

Author: Benjamin Jean-Marie Tremblay

Date: 11 July 2021

## Software info

All code was run using R version 4.1.0 (aarch64-apple-darwin20.4.0, 64-bit).
The follow packages were also used:

### From CRAN

- `ggplot2` v3.3.3
- `sf` v0.9.8
- `rnaturalearth` v0.1.0
- `rnaturalearthdata` v0.1.0
- `dplyr` v1.0.6
- `viridis` v0.6.1
- `cowplot` v1.1.1
- `circlize` v0.4.12
- `tidyr` v1.1.3
- `readr` v1.4.0
- `reshape2` v1.4.4

### From Bioconductor (version 3.13)

- `Biostrings` v2.60.1
- `GenomicRanges` v1.44.0
- `ComplexHeatmap` v2.8.0
- `Rsamtools` v2.8.0
- `GenomicFeatures` v1.44.0

## Files


- `1.map.R`: Generate the map and timescale plots
- `2.circos.R`: Generate the percent ID circos plots and a TSV of sample
  percent IDs
- `3.toxin.R`: Generate the toxin plot and TSV of toxin SNPs
- `4.compareToxinVariants.R`: Compare the toxin MSA with Mike's toxin MSA
- `5.circosLegend.R`: Make a nicer legend for the circos percent ID plots
- `maskRef.R`: Example script to mask the reference bases in sequences
- `plotCoverage.R`: A small program to generate coverage circos from BAM files
- `plotCoverage-plasmid.R`, `plotCoverage-chromosome.R`: Modifications of the
  program that are to be used with the chromosome and plasmid BAM files so
  that they end up looking similar to the ones produced by `2.circos.R`
- `6.plotCoverage-plasmid.sh`, `plotCoverage-chromosome.sh`: Bash scripts to
  execute the above programs. The plasmid one can be executed from this
  directory, but the chromosome one was run from the `dox3` machine.

