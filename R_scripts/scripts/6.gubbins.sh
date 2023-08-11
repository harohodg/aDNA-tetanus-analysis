# Benjamin Jean-Marie Tremblay

set -ex

run_gubbins \
  --starting_tree data/core-snippy-RAxML_tree.tre \
  --verbose \
  --outgroup cochlearium \
  --prefix results/gubbins/core-snippy-genome \
  --threads 8 \
  --filter_percentage 90 \
  data/core-snippy-genome.aln

