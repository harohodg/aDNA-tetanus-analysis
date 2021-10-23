#!/bin/bash

SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SINGULARITYDIR="${SCRIPTDIR}"/../tools/
REFGENOME="${SCRIPTDIR}"/../data/ref_genome/NC_004565.1_plasmid_genomic_revComplement.fna
REFTENTSEQS="${SCRIPTDIR}"/../data/ref_genome/Reference_TeNT_sequences.fasta

mkdir -p "${SCRIPTDIR}"/../data/TeNT
cd "${SCRIPTDIR}"/../data/TeNT

if [[ $# -eq 0 ]]
then
	echo
	echo "Usage: variants.sh --variants --alignments --trees"
	echo "    Input options:"
	echo "    -v, --variants       Start analysis at variant calling."
	echo "    -a, --alignments     Start analysis at alignment manipulation."
	echo "    -t, --threads        Number of threads to use."
	echo "                         Default: all threads"
	echo "    -h, --help           Print this help message."
	echo
	exit 1
fi

POSITIONAL=()
while [[ $# -gt 0 ]];
do
	key="$1"
	case $key in
		-v|--variants)
			VARIANTS="TRUE"
			shift
		;;
		-a|--alignments)
			ALIGNMENTS="TRUE"
			shift
		;;
		-t|--threads)
			THREADS="$2"
			shift
			shift
		;;
		-h|--h|-help|--help)
			echo
			echo "Usage: variants.sh --variants --alignments --trees"
			echo "    Input options:"
			echo "    -v, --variants       Start analysis at variant calling."
			echo "    -a, --alignments     Start analysis at alignment manipulation."
			echo "    -t, --threads        Number of threads to use. Default: all threads."
			echo "                         Note that the parallelization is not particularly efficient, so be patient!"
			echo "                         Default: all threads"
			echo "    -h, --help           Print this help message."
			echo
			exit 1
		;;
		*)
			# Takes any other argument, then throws it out
			POSITIONAL+=("$1")
			shift
		;;
	esac
done

set -- "${POSITIONAL[@]}"

# If it's given an option it doesn't recognize, just stop and fail
if [[ -n $1 ]]; then
	echo
	echo "    Error: Unknown argument or option: $1 $2"
	echo "    Use variants.sh -h for a list of available options."
	echo
	exit 1
fi

# If num cpus was not specified, it's the full system. If it was specified, make sure it is an integer.
if [ -z ${THREADS+x} ]
then
	THREADS=$(nproc)
elif ! [[ $THREADS =~ ^-?[0-9]+$ ]]
then
	echo
	echo "    Error: invalid option for -t/--threads: $THREADS"
	echo
	exit 1
fi

if [ "${VARIANTS}" == "TRUE" ];
then
	mkdir -p misc
	if [ ! -f ./misc/NC_004565.1_revComplement.fasta.bwt ];
	then
		ln -srf "${REFGENOME}" ./misc/NC_004565.1_revComplement.fasta

		singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools faidx ./misc/NC_004565.1_revComplement.fasta
		# I just used BLAST to figures out gene coordinates as their reverse
		# complemented ones.
		singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools faidx ./misc/NC_004565.1_revComplement.fasta NC_004565.1:1496-5443 > ./misc/tmp.fasta

		mv ./misc/tmp.fasta ./misc/NC_004565.1_revComplement_tetX.fasta
		singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools faidx ./misc/NC_004565.1_revComplement_tetX.fasta
		singularity exec "${SINGULARITYDIR}"/bwa_0.7.17--h5bf99c6_8.sif bwa index ./misc/NC_004565.1_revComplement.fasta
	fi

	# Iterate over a bunch of parameters.
	# The first sets the minimum posterior probability of the variant.
	for POSTERIOR in 95;
	do
	POSTERIOR="0."${POSTERIOR}"" # The posterior is a percentage. ;)
		# The second is the minimum fraction of "good" bases at a position.
		for MINGOODFRAC in 75;
		do
		MINGOODFRAC="0."${MINGOODFRAC}"" # again, as percentage
			# The third is the minimum pileup quality for a position to be considered.
			for MINPILEUP in 35;
			do
				for SET in merged singles;
				do
					mkdir -p "${SET}"
					cd "${SET}"

					find ../../read_alignments/plasmid/"${SET}" -iname "*.bam" | while read BAM;
					do
						# The Octopus variant caller is sensitive to read groups indicated in the alignment.
						# I'm extracting the reads and realigning them with manual RG tags.
						BASENAME=$(basename "${BAM}" | cut -d "_" -f 1)
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools view -bS -@ "${THREADS}" -F 4 -o tmp.bam "${BAM}"
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools sort -@ "${THREADS}" -o tmp.sort.bam tmp.bam
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools index tmp.sort.bam
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools fastq tmp.sort.bam > "${BASENAME}".fastq
						rm tmp.bam tmp.sort.bam tmp.sort.bam.bai
						singularity exec "${SINGULARITYDIR}"/bwa_0.7.17--h5bf99c6_8.sif bwa mem -t "${THREADS}" -R "@RG\tID:1\tSM:"${BASENAME}"\tLB:octo\tPU:illumina" ../misc/NC_004565.1_revComplement.fasta "${BASENAME}".fastq > "${BASENAME}".sam
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools view -@ "${THREADS}" -bh -o "${BASENAME}".unsorted.bam "${BASENAME}".sam
						rm "${BASENAME}".sam "${BASENAME}".fastq
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools sort -@ "${THREADS}" -o "${BASENAME}".sorted.bam "${BASENAME}".unsorted.bam
						rm "${BASENAME}".unsorted.bam
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools index "${BASENAME}".sorted.bam
						singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif samtools depth -a "${BASENAME}".sorted.bam > "${BASENAME}".depth

						singularity exec "${SINGULARITYDIR}"/octopus_0.7.4.sif \
							octopus \
								--threads "${THREADS}" \
								--reference ../misc/NC_004565.1_revComplement.fasta \
								--reads "${BASENAME}".sorted.bam \
								--allow-octopus-duplicates \
								--regions NC_004565.1:1496-5443 \
								--output "${BASENAME}".octopus.vcf \
								--bamout "${BASENAME}".octopus.bam \
								--mask-low-quality-tails 5 \
								--min-mapping-quality 10 \
								--min-pileup-base-quality "${MINPILEUP}" \
								--min-good-base-fraction "${MINGOODFRAC}" \
								--min-variant-posterior "${POSTERIOR}" \
								--organism-ploidy 1

							rm "${BASENAME}".sorted.bam "${BASENAME}".sorted.bam.bai

							singularity exec "${SINGULARITYDIR}"/bcftools_1.12--h45bccc9_1.sif bcftools view --output-type z --output "${BASENAME}".octopus.vcf.gz "${BASENAME}".octopus.vcf
							singularity exec "${SINGULARITYDIR}"/samtools_1.12--h9aed4be_1.sif htsfile "${BASENAME}".octopus.vcf.gz
							singularity exec "${SINGULARITYDIR}"/bcftools_1.12--h45bccc9_1.sif bcftools index "${BASENAME}".octopus.vcf.gz
							# Get consensus calls. Note that reference positions are masked in lowercase.
							singularity exec "${SINGULARITYDIR}"/bcftools_1.12--h45bccc9_1.sif bcftools consensus --fasta-ref ../misc/NC_004565.1_revComplement_tetX.fasta --output "${BASENAME}".consensus.fasta "${BASENAME}".octopus.vcf.gz
							# Normally, you could do simple VCF filtering like: --include 'MIN(FMT/DP>5)'.
							# In this case, we have to do something a little fancier to report the desired results.
							# Instead of filtering the VCF file before consensus calling, I do some post-filtering with python to ensure
							# no-coverage positions are masked properly.
							# This can be changed to allow for many different filtering options.
							python "${SCRIPTDIR}"/mask_low_coverage_positions.py --fastaFile "${BASENAME}".consensus.fasta --depthFile "${BASENAME}".depth --region 1496-5443 --outFile "${BASENAME}".consensus.masked.fasta --noCoverageChar - --lowCoverageChar lower
					done
					cd ../
				done

				mkdir -p results/summary results/archive results/variants

				# Make a nice table.
				find merged singles -iname '*vcf' | while read VARIANTFILE;
				do
					NUMVAR=$(grep -v "#" "${VARIANTFILE}" | wc -l)
					BASENAME=$(basename "${VARIANTFILE}" | awk -F '.octopus' '{print $1}')
					printf "$BASENAME\t$POSTERIOR\t$MINGOODFRAC\t$MINPILEUP\t$NUMVAR\n"
				done > summary.minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".tsv

				# Merge the VCF files
				singularity exec "${SINGULARITYDIR}"/bcftools_1.12--h45bccc9_1.sif bcftools merge --output-type v --output ./merged_variants.minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".vcf --threads "${THREADS}" ./merged/*.vcf.gz ./singles/*.vcf.gz
				find merged singles -iname '*.vcf' | while read i; do b=$(basename $i | awk -F '.octopus' '{print $1}'); grep -v '^#' $i | while read v; do printf "$b\t$v\n"; done; done > all_variants.long.tsv
				find merged -iname '*.vcf' | while read i; do b=$(basename $i | awk -F '.octopus' '{print $1}'); grep -v '^#' $i | while read v; do printf "$b\t$v\n"; done; done > merged_variants.long.tsv


				# Summarize
				env GZIP=-9 tar czf minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".tar.gz ./merged ./singles
				mv minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".tar.gz ./results/archive
				mv summary.minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".tsv ./results/summary
				mv ./merged_variants.minPosterior_"${POSTERIOR}".minGoodFrac_"${MINGOODFRAC}".minPileupQual_"${MINPILEUP}".vcf all_variants.long.tsv merged_variants.long.tsv ./results/variants
			done
		done
	done
fi

if [ "${ALIGNMENTS}" == "TRUE" ];
then
	# First, align nucleotide sequences to the reference
	# sequence using mafft --addfragments.
	find merged -iname '*consensus.masked.fasta' | while read FASTAFILE;
	do
		BASENAME=$(basename "${FASTAFILE}" | awk -F '.fasta' '{print $1}')
		SEQ=$(grep -v "^>" "${FASTAFILE}")
		printf ">$BASENAME\n$SEQ\n"
	done > tmp.fasta

	# Exclude sequences with less than 1 position.
	cat tmp.fasta | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | while read LINE;
	do
		if [[ "${LINE}" == *">"* ]];
		then
			HEADER="${LINE}"
		else
			SEQ="${LINE}"
			SEQLEN=$(echo "${SEQ}" | tr -d '\n' | sed 's/-//g' | wc -m)
			if [ "${SEQLEN}" -ge 1 ];
			then
				printf "$HEADER\n$SEQ\n"
			fi
		fi
	done > consensus_variants.fasta
	rm tmp.fasta

	# Make a reference alignment from reference sequences.
	singularity exec "${SINGULARITYDIR}"/mafft_7.480--h779adbc_0.sif mafft --auto --thread "${THREADS}" "${REFTENTSEQS}" > reference_tents.linsi.fasta

	# Align consensus variant fragments to reference alignment.
	# Note that the default behaviour allows for insertions in the variants, and because of this some sequences end up out of frame.
	singularity exec "${SINGULARITYDIR}"/mafft_7.480--h779adbc_0.sif mafft --auto --thread "${THREADS}" --addfragments consensus_variants.fasta reference_tents.linsi.fasta > consensus_variants.E88_aligned.withInsertions.fasta
	singularity exec "${SINGULARITYDIR}"/mafft_7.480--h779adbc_0.sif mafft --auto --keeplength --mapout --thread "${THREADS}" --addfragments consensus_variants.fasta reference_tents.linsi.fasta > consensus_variants.E88_aligned.noInsertions.fasta

	# Second, make protein alignments.
	# I use a custom script here, because looking at the nucleotide alignment above,
	# it is basically an in-frame codon alignment as it is. So, instead of translating
	# to protein and trying to realign the sequences, I do an in-frame translation
	# of each codon in the variant sequences, and highlight synonymous and non-
	# synonymous substitutions.
	python "${SCRIPTDIR}"/translate_in_frame.py --fastaFile consensus_variants.E88_aligned.noInsertions.fasta --outputFasta consensus_variants.E88_aligned.noInsertions.AA.fasta --outputTSV consensus_variants.E88_aligned.noInsertions.AA.tsv

	# Summarize
	mkdir -p results/fasta
	mv *.fasta consensus_variants.fasta.map consensus_variants.E88_aligned.noInsertions.AA.tsv ./results/fasta
fi
