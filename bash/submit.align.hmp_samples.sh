#!/bin/bash

#Define base files and folders
useSamples=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/parameter_files/sample_list.combined;
script=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/shell/align.hmp_samples.sh;

#Submit for parallel processing
parallel --no-notice -j20 --tag -a $useSamples $script;
