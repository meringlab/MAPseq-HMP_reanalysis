#!/bin/bash
################################################################################
#Global reference-based OTU demarcation using UCLUST
#=> at 97% (level "R")
#=> mapping against a (UCLUST-preclustered) database of representative sequences
#
#Export:
#=> OTU files
#=> sequence-to-OTU mapping
#=> OTU-to-sequence mapping
#=> OTU datas
#=> OTU tables
#
#
#2015-11-20
#sebastian.schmidt@imls.uzh.ch
################################################################################

#Set path to uclust binary
UCLUST=~jfmrod/usr/bin/uclust;
SORTABUND=/local/erisdb/sebastian/toolbox_16S/2-perl/sort_fasta_by_abundance.pl;
DEREPLICATE=/local/erisdb/sebastian/toolbox_16S/2-perl/dereplicate_fasta.pl;

#Set basic parameters to values used in QIIME HMP SOP
threshold=0.97                    #clustering threshold

#Set parameters to QIIME defaults
#=> according to https://qiime.wordpress.com/2010/12/17/new-default-parameters-for-uclust-otu-pickers/
maxaccepts=20                     #maximum hits to look through
maxrejects=500                    #maximum potential hits to reject before matching
wordlength=12                     #length of kmer for searching

#Set basic directory and file names
export folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
export folderData=$folderBase/1-data/global;
export folderOutput=$folderBase/1-data/global;
export folderRefDB=/local/erisdb/sebastian/toolbox_16S/reference_databases/full-ssu-db-ref;
export fileDB=$folderRefDB/full-ssu-db3.clustering.bacteria.uc.R.rep.fasta;

#Turn "U" characters in databaset to "T", because UCLUST doesn't otherwise understand what it's supposed to do
cat $fileDB | perl -pi -e 's/U/T/g' > tmp.db;
mv tmp.db $fileDB;

#Iterate over subregions
for v in v13 v35; do
  export V=$v;
  #Set current file names
  
  export fileSampleMapping=$folderData/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileUnaligned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.fasta;
  export fileUnalignedPruned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.pruned.fasta;
  export fileDerepFASTA=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.derep.fasta;
  export fileDerepCluster=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.derep.cluster;
  export fileUnalignedSortedRC=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.sorted.rev_comp.fasta;
  export fileSampleMapping=$folderOutput/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileTaxonomy=$folderOutput/$v.nonchimeric.filtered.full_db.taxonomy;
  export fileSeqAccHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_id.hash;
  export fileClosedUC=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.uc;
  export fileClosedRefOTU=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.otu;
  export fileClosedRefMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.r.otu;
  export fileClosedRefSeqMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.r.seq_mapping.otu;
  export fileClosedRefOTUTable=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.r.otu_table.tsv;
  export fileClosedRefOTUData=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.closed_ref.r.otu_data.tsv;
  
  #For V13 sequences, chop off leading ~70bp to obtain a sequence that UCLUST can match to the database
  if [ "$v" = "v13" ]; then
    cat $fileUnaligned | perl -ane 'if (/^>/) {print} else {chomp; $seq=$_; $l=length $seq; $start=$l-400; if($start<=0){print "$seq\n"} else{$keep=substr $seq,$start; print "$keep\n"}}' > $fileUnalignedPruned;
  fi
  
  #Dereplicate FASTA file
  $DEREPLICATE $fileUnalignedPruned $fileDerepFASTA $fileDerepCluster;
  
  #Run UCLUST on sorted sequences
  $UCLUST --input $fileDerepFASTA --id 0.97 --lib $fileDB --libonly --rev --maxaccepts $maxaccepts --maxrejects $maxrejects --w $wordlength --uc $fileClosedUC;
  
  
done









