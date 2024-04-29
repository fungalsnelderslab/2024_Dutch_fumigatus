These scripts are used to align the reads, do the read mapping and variant calling/filtering

The aligner.sh file is in common for all samples. For the filtering, variant calling, and analysis the scripts are split into two sections. The  `dutch_filtering.sh` and `dutch_analysis.sh` files are used for analysis within the Dutch population (PCA, Fst, phylogenetic networks). The `global_filtering.sh` and `global_analysis.sh` scripts are used to produce the data needed for all analyses on the wider dataset.
