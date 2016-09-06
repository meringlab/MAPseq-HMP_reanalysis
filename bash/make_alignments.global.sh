#!/bin/bash
################################################################################
#Make global alignments from per-sample alignments
#=> exclude samples with <1000 sequences (after filtering)
#=> concatenate per-sample alignments
#=> de-replicate sequences, accepting ~1% mismatches against abundant "seeds"
#
#Export (R-readable output files):
#=> global alignments (de-replicated)
#=> global unaligned files (de-replicated)
#=> global sequence-to-sample mapping
#
#
#2015-08-11
#sebastian.schmidt@imls.uzh.ch
################################################################################

#Set path to "denoiser" binary
DENOISER=/home/jfmrod/gaia/usr/bin/denoiser;

#Set basic parameters
threads=80                        #number of threads to use in denoising
maxMismatches=5                   #maximum accepted mismatches of aligned sequence to abundant "seed"
minSampleSize=1000                #minimum number of (filtered) sequences per sample
export minSeqFreq=3               #minimum global abundance per sequence after denoising
export minSamplePenetrance=2      #minimum number of unique samples that a sequence must be observed in, after denoising

#Set basic directory and file names
folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
folderData=$folderBase/1-data/per_sample;
folderOutput=$folderBase/1-data/global;

#Iterate over subregions
for v in v13 v35 v69; do
  #Get current sample lists
  useSamples=$folderBase/parameter_files/sample_list.$v;
  fileSamples=$folderBase/parameter_files/sample_list.$v.tmp;
  fileSamplesFiltered=$folderBase/parameter_files/sample_list.$v.filtered;
  cut -f 1 $useSamples > $fileSamples;
  
  #Remove all samples that contain too few sequences
  cat $fileSamples | while read f; do
    size=$(gunzip -c $folderData/$f/$v.fasta.gz | grep -c "^>");
    if [ "$size" -gt "$minSampleSize" ]; then echo $f; fi
  done > $fileSamplesFiltered;
  
  #Generate global temporary alignment
  export fileAlignmentRAW=$folderOutput/$v.nonchimeric.filtered.aligned.global.raw.fasta;
  echo -n "" > $fileAlignmentRAW;
  cat $fileSamplesFiltered | while read f; do gunzip -c $folderData/$f/$f.$v.nonchimeric.filtered.aligned.sto.gz | perl -ane '($acc, $seq)=split / /; print ">$acc\n$seq" unless $acc =~ /\/\//' >> $fileAlignmentRAW; done
  
  #Generate global unaligned file
  export fileUnaliged=$folderOutput/$v.nonchimeric.filtered.unaligned.global.fasta;
  echo -n "" > $fileUnaliged;
  cat $fileSamplesFiltered | while read f; do gunzip -c $folderData/$f/$f.$v.nonchimeric.filtered.fasta.gz >> $fileUnaliged; done
  
  #Run denoiser to de-replicate alignment
  export fileDenoised=$folderOutput/$v.nonchimeric.filtered.aligned.global.raw.denoised.cluster;
  $DENOISER -aligned true -maxmiss $maxMismatches -nthreads $threads -of tmp $fileAlignmentRAW > $fileDenoised;
  
  #Filter data
  export fileDenoisedFASTA=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.fasta;
  export fileDenoisedUnaligned=$folderOutput/$v.nonchimeric.filtered.unaligned.global.denoised.fasta;
  export fileSeqOrder=$folderOutput/$v.nonchimeric.filtered.denoised.seq;
  export fileSeqAccHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_id.hash;
  perl -e '
    use Storable;
    
    #Read denoised clustering & filter sequences by minimum abundance criterion
    open(CLUSTER, $ENV{fileDenoised}) or die;
    while (<CLUSTER>) {
      next if /^#/;
      chomp; ($seed, $size, @children) = split /\t/;
      #print "Removed $seed\t$size\n" if ($size <= ($ENV{minSeqFreq} - 1));
      next if ($size <= ($ENV{minSeqFreq} - 1));
      push @seq_order => $seed;
      foreach $child (@children) {$Seed{$seed}{Children}{$child}=1; $Child{$child}=$seed; $Child{$seed}=$seed; push @seq_order => $child}
    }
    close CLUSTER;
    
    #Read FASTA file and store seed sequences
    open(FASTA, $ENV{fileAlignmentRAW}) or die;
    while ($acc = <FASTA>) {
      chomp $acc;
      if ($acc =~ /^>/) {$acc =~ s/^>//; if (exists $Seed{$acc}) {$seq=<FASTA>; chomp $seq; $Seed{$acc}{Seq}=$seq}}
    }
    close FASTA;
    
    #Read unaligned file
    open(UNA, $ENV{fileUnaliged}) or die;
    while ($acc = <UNA>) {
      chomp $acc;
      if ($acc =~ /^>/) {
        $acc =~ s/^>//;
        $seq = <UNA>; chomp $seq;
        $Unaligned{$acc} = $seq;
      }
    }
    close UNA;
    
    #Export sequences & list of retained sequences
    open(OUT, ">$ENV{fileDenoisedFASTA}") or die;
    open(UNAOUT, ">$ENV{fileDenoisedUnaligned}") or die;
    open(SEQ, ">$ENV{fileSeqOrder}") or die;
    $id = 1;
    foreach $acc (@seq_order) {print OUT ">$acc\n$Seed{$Child{$acc}}{Seq}\n"; print UNAOUT ">$acc\n$Unaligned{$acc}\n"; print SEQ "$acc\n"; $ID{$id}=$acc; $Acc{$acc}=$id; $id++}
    close OUT;
    close SEQ;
    
    #Store Seq_Acc and Seq_ID hashes
    store \%Acc, $ENV{fileSeqAccHash}; store \%ID, $ENV{fileSeqIDHash};
  ';
  
  #Generate global, filtered sequence->sample mapping file
  export fileMappingTMP=$folderOutput/$v.mapping.tmp;
  export fileMapping=$folderOutput/$v.nonchimeric.filtered.denoised.seq2sample.map;
  echo -n "" > $fileMappingTMP;
  cat $fileSamplesFiltered | while read f; do export F=$f; gunzip -c $folderData/$f/$f.$v.nonchimeric.filtered.aligned.sto.gz | perl -ane '($acc, $seq)=split / /; print "$acc\t$ENV{F}\n" unless $acc =~ /\/\//' >> $fileMappingTMP; done
  perl -e '
    use Storable;
    %Acc = %{retrieve($ENV{fileSeqAccHash})};
    
    #Read temporary mapping and export mapping for allowed sequences on the fly
    open(MAP, $ENV{fileMappingTMP}) or die;
    open(OUT, ">$ENV{fileMapping}") or die;
    while ($line=<MAP>) {($acc, $smpl)=split /\t/ => $line; print OUT $line if exists $Acc{$acc}}
  ';
  rm $fileMappingTMP;
done









