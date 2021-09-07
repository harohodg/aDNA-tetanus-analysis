# Benjamin Jean-Marie Tremblay
# 2021-07-11

scripts/./plotCoverage-plasmid.R \
  --bam=data/bam-plasmid/ \
  --gff=data/E88_plasmid.gff3 \
  --bin-size=300 \
  --hl-features=NC_004565.1:68640..72587 \
  --hl-colours=red \
  --track-order=results/plotCoverage-plasmid_Order.tsv \
  --out=figures/CircosCoverage_plasmid_SAM.pdf

scripts/./plotCoverage-plasmid.R \
  --bam=data/bam-plasmid/ \
  --gff=data/E88_plasmid.gff3 \
  --bin-size=300 \
  --plot-sum=true \
  --hl-features=NC_004565.1:68640..72587 \
  --hl-colours=red \
  --track-order=results/plotCoverage-plasmid_Order.tsv \
  --out=figures/CircosCoverage_plasmid_SAM_totalCoverage.pdf
