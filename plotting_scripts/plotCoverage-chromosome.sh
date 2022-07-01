# Benjamin Jean-Marie Tremblay
# 2021-07-11

./plotCoverage-chromosome.R \
  --bam=/disk3/acdoxey/aDNA-tetani/fastp_trimmed_E88_bowtie2_mapping_pooled_files/chromosome/merged/ \
  --gff=/disk3/acdoxey/aDNA-tetani/fastp_trimmed_E88_bowtie2_mapping_pooled_files/E88_chromosome.gff3 \
  --bin-size=11250 \
  --track-order=plotCoverage-chromosome_Order.tsv \
  --R-lib=R_lib \
  --out=CircosCoverage_chromosome_SAM.pdf

./plotCoverage-chromosome.R \
  --bam=/disk3/acdoxey/aDNA-tetani/fastp_trimmed_E88_bowtie2_mapping_pooled_files/chromosome/merged/ \
  --gff=/disk3/acdoxey/aDNA-tetani/fastp_trimmed_E88_bowtie2_mapping_pooled_files/E88_chromosome.gff3 \
  --bin-size=11250 \
  --plot-sum=true \
  --track-order=plotCoverage-chromosome_Order.tsv \
  --R-lib=R_lib \
  --out=CircosCoverage_chromosome_SAM_totalCoverage.pdf
