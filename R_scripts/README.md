# R analysis and plotting scripts

## `scripts/`

Analysis and plotting scripts for various figures. These were run using R version 4. All packages used can be installed from CRAN and Bioconductor. The Gubbins script was run using Gubbins version 3.3.0.

## `data/`

Some input data (anything <1 MB) required for the scripts to work can be found in the `data/` folder. The remaining data can be obtained as the output of previous analyses (and deposited in the `data/` folder to make them available for the R scripts). These include:

* `core.full.cleaned.aln`: Alignment file from previous analyses
* `core-snippy-genome.aln`: Alignment file from previous analyses
* `bam-plasmid/*.bam`: BAMs of reads aligned to the E88 plasmid from previous analyses
* `E88_chromosome.gff3`: E88 chromosome sequence annotation data from NCBI
* `E88_plasmid.gff3`: E88 plasmid sequence annotation data from NCBI

See the NCBI page for E88 [here](https://www.ncbi.nlm.nih.gov/genome/1098?genome_assembly_id=300529).

