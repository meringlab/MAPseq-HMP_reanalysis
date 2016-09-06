#!/bin/bash
################################################################################
#Align all nonchimeric sequences, using Infernal
#=> cut at different positions, depending on the current subregion of choice
#=> filter out sequences which do not align to the bacterial model, align too shortly or with too many gaps
#
#Export:
#=> trimmed alignment
#
#References:
#Infernal:
#Nawrocki et al, Bioinformatics, 2009
#Nawrocki et al, Bioinformatics, 2013
#ssu-align (source for the 16S rRNA model)
#Nawrocki, PhD Thesis, Washington U St Louis, 2009
#
#Online documentation for Infernal:
#http://infernal.janelia.org
#
#2015-08-06
#sebastian.schmidt@imls.uzh.ch
################################################################################

f=$1;

#Set path to Infernal binary and bacterial 16S model
INFERNAL=/home/jfmrod/erisdb/usr/src/infernal-1.1rc4/src/cmalign;
MODEL=/local/atlas/sebastian/toolbox_16S/infernal_models/bacteria-0p1.1.cm

#Define flanking positions
export v13start=30;
export v13end=539;
export v35start=410;
export v35end=933;
export v69start=1039;
export v69end=1530;

#Set basic parameters
export maxGaps=100                               #Accept ~20% gaps
export maxInserts=100                            #Accept ~20% inserts (unaligned bases)
export maxLeadingGaps=50                         #Accept maximum 50 leading gaps

#Set basic directory and file names
folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
folderData=$folderBase/1-data/per_sample;
folderSample=$folderData/$f;
folderOutput=$folderSample/alignment;
mkdir $folderOutput;


#Iterate through subregions
for v in v13 v35 v69; do
  export fileNC=$folderSample/$f.$v.nonchimeric.fasta;
  if [ -e $fileNC.gz ]; then
    export V=$v;
     
    #Gunzip
    gunzip -f $fileNC.gz;
    
    #Get current file names
    fileAlignedTMP=$folderOutput/tmp.align.$v;
    export fileTMPSTO=$folderOutput/tmp.align.$v.sto;
    export fileAlignLog=$folderOutput/tmp.align.$v.log;
    #Output file names
    export fileAlignScores=$folderOutput/$f.$v.scores;
    export fileAlignedSTO=$folderSample/$f.$v.nonchimeric.filtered.aligned.sto;
    export fileUnalignedFiltered=$folderSample/$f.$v.nonchimeric.filtered.fasta;
    export fileFilteredSeq=$folderSample/$f.$v.nonchimeric.filtered.seq;
    export SeqAccHash=$folderSample/$f.$v.nonchimeric.filtered.seq_acc.hash;
    export SeqIDHash=$folderSample/$f.$v.nonchimeric.filtered.seq_id.hash;
    fileAlignReport=$folderOutput/$f.$v.report;
    
    #Run Infernal
    $INFERNAL --sub --notrunc --cpu 40 --outformat Pfam -g -o $fileAlignedTMP $MODEL $fileNC | egrep -v '^(#|$)' > $fileAlignLog;
    cat $fileAlignedTMP | egrep -v '^(# |$|#=|//$)' > $fileTMPSTO;
    rm $fileAlignedTMP;
    
    #Process output
    perl -e '
      use Storable;
      
      #Preallocate
      my %Sequences;
      my %Filtered;
      my %IDs;
      
      #Process scores
      open(RAW, $ENV{fileAlignLog});
      open(SCORES, ">$ENV{fileAlignScores}");
      while (<RAW>) {
        chomp;
        #Edit whitspace (spaces to tabs)
        s/ +/\t/g;
        @line = split /\t/;
        
        #Require positive alignment score against bacterial model
        if ($line[7] > 0 & $line[12] > 0) {
          print SCORES "$line[2]\tbacteria\t$line[7]\t$line[12]\n";
          $Sequences{$line[2]}{Score} = $line[7];
        }
        else {$Sequences{$line[2]}{Score} = $line[7]; print SCORES "$line[2]\tnone\t0\t0\n"}
      }
      close RAW; close SCORES;
      
      #Process sequences
      open(SEQ, $ENV{fileTMPSTO});
      open(ALIGNED, ">$ENV{fileAlignedSTO}");
      open(SEQID, ">$ENV{fileFilteredSeq}");
      $idx=1;
      LINE:
      while (<SEQ>) {
        chomp;
        #Edit whitspace (spaces to tabs)
        s/ +/\t/g;
        #Get current sequence
        ($acc, $alignment) = split /\t/;
        $tmp = $alignment;
        
        #Skip sequence if score is too low
        if ($Sequences{$acc}{Score} < 0) {print "Sequence has negative alignment score:\t$acc\t$alignment\n"; next LINE}
        
        #Skip sequences with ³20% unaligned bases (i.e., >=100 insertions)
        if (($tmp =~ tr/acgun//) >= $ENV{maxInserts}) {print "Sequence has to many inserts:\t$acc\t$alignment\n"; next LINE}
        
        #Remove "points" (gaps in sequence)
        $alignment =~ s/\.//g;
        #Remove insertions
        $alignment =~ s/[a-z]//g;
        
        #Prune alignment
        my $cutStart; my $cutEnd;
        if ($ENV{V} eq "v13") {$cutStart=$ENV{v13start}; $cutEnd=$ENV{v13end}}
        elsif ($ENV{V} eq "v35") {$cutStart=$ENV{v35start}; $cutEnd=$ENV{v35end}}
        elsif ($ENV{V} eq "v69") {$cutStart=$ENV{v69start}; $cutEnd=$ENV{v69end}}
        
        $start = $cutStart - 1; $length = $cutEnd - $cutStart;
        $pruned = substr $alignment, $start, $length;
        
        #Skip if sequence is too gappy
        $tmp = $pruned;
        if (($tmp =~ tr/-//) >= $ENV{maxGaps}) {print "Sequence too gappy:\t$acc\t$pruned\n"; next LINE}
        
        #Skip if sequence has two (or more) consecutive Ns
        if ($pruned =~ /N-*N/) {print "Sequence has too many Ns:\t$acc\t$pruned\n"; next LINE}
        
        #Skip if sequences starts with too many gaps (10 or more)
        #if ($pruned =~ /^-{$ENV{maxLeadingGaps},}/) {print "Sequence has too many leading gaps:\t$acc\t$pruned\n"; next LINE}
        
        #Export
        print ALIGNED "$acc $pruned\n";
        print SEQID "$idx\t$acc\n";
        
        $Filtered{$acc} = $idx;
        $IDs{$idx} = $acc;
        
        $idx++;
      }
      print ALIGNED "//\n";
      close SEQID; close ALIGNED; close SEQ;
      
      
      #Export unaligned sequences
      open(FASTA, "$ENV{fileNC}");
      open(UNALIGNED, ">$ENV{fileUnalignedFiltered}");
      LINE:
      while ($line = <FASTA>) {
       chomp $line;
       if ($line =~ /^>/) {
         $line =~ s/^>//;
         $seqline = <FASTA>;
         
         if ($Filtered{$line} > 0) {print UNALIGNED ">$line\n$seqline"}
        }
      }
      close FASTA; close UNALIGNED;
      
      
      #Export sequence ID and Acc hashes
      #V13
      my %Acc = (); my %ID = ();
      foreach $id (sort {$a <=> $b} keys %IDs) {$ID{$id} = $IDs{$id}; $Acc{$IDs{$id}} = $id}
      store \%Acc, $ENV{SeqAccHash}; store \%ID, $ENV{SeqIDHash};
    ' > $fileAlignReport;
    
    #Tidy up
    gzip -f $fileTMPSTO;
    gzip -f $fileAlignLog;
    gzip -f $fileAlignScores;
    gzip -f $fileAlignedSTO;
    gzip -f $fileUnalignedFiltered;
    gzip -f $fileNC;
  fi
done

