# MAPseq-HMP_reanalysis
This repository contains code used to reanalyse Human Microbiome Project data to benchmark MAPseq.

Raw data from the HMP is pre-processed and matched to samples and sequencing subregions in the script [prepare.raw_data.pl](perl/prepare.raw_data.pl). This generates one folder per sample for subsequent sample-wise (parallel) processing. It also generates global mapping tables.

The script [submit.filter_chimeras.sh](bash/submit.filter_chimeras.sh) then calls the script [filter.chimeric_reads.sh](bash/filter.chimeric_reads.sh) for each sample in parallel to remove chimeras.

The script [submit.align.hmp_samples.sh](bash/submit.align.hmp_samples.sh) calls the script [align.hmp_samples.sh](bash/align.hmp_sample.sh) for each sample in parallel to run INFERNAL. Afterwards, the script [make_alignments.global.sh](bash/make_alignments.global.sh) stitches together individual per-sample alignments and de-replicates and de-noises them a bit. This script also defines the "consensus" lists of sequences to be used downstream by clustering methods – keeping only those sequences which are non-chimeric and align satisfactorily to the target subregion.

The scripts "make_otus.*.sh" then call the different mapping/clustering tools to generate OTU sets and translate them into R-readable formats (OTU tables etc.).

Finally, the script [hmp.benchmark.R](R/hmp.benchmark.R) contains all the code for R-based analyses, as detailed in the manuscript.

(Raw) results are available in the folder [results](results).
