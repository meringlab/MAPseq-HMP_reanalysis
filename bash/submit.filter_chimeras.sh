#!/bin/bash

#Define base files and folders
useSamples=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/parameter_files/sample_list.combined;
script=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/shell/filter.chimeric_reads.sh;

#Submit for parallel processing
parallel --no-notice -j40 --tag -a $useSamples $script;
