#!/bin/bash
################################################################################
#Filter for chimeric sequences
#=> run UCHIME in reference-based mode
#=> against our in-house reference database of nonchimeric sequences
#=> run UCHIME in de_novo mode
#=> combine outputs of both modes
#
#Export:
#=> nonchimeric sequence FASTA files
#
#2015-08-06
#sebastian.schmidt@imls.uzh.ch
################################################################################

#Get current sample name (as passed from GNU Parallel)
f=$1;

#Get binaries
USEARCH=/local/atlas/sebastian/bin/usearch;

#Set parameters
export minLength=400;

#Get current files and folders
folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
folderData=$folderBase/1-data/per_sample;
folderSample=$folderData/$f;
ChimeraReference=/mnt/mnemo2/sebastian/projects/1-otu_testing/joao.reference.complete.uc.new.fasta;
ChimeraReferenceUDB=/mnt/mnemo2/sebastian/projects/1-otu_testing/joao.reference.complete.uc.new.udb;    #pre-built index

#Set output folder and tidy up
folderFilterChimera=$folderSample/filter.chimera
mkdir $folderFilterChimera;
rm -f $folderFilterChimera/*;

#Iterate through subregions
for v in v13 v35 v69; do
  export fileRAW=$folderSample/$v.fasta;
  #Gunzip, just to be sure
  gunzip $fileRAW.gz;
  if [ -e $fileRAW ]; then
    #Set current file names
    export fileOutUCRef=$folderFilterChimera/$f.$v.uchime_ref.uc;
    export fileOutChimerasRef=$folderFilterChimera/$f.$v.uchime_ref.chimeras.fasta;
    export fileOutNonchimerasRef=$folderFilterChimera/$f.$v.uchime_ref.nonchimeric.fasta;
    export fileDerepOut=$folderFilterChimera/$f.$v.derep.uc;
    export fileDerepFASTA=$folderFilterChimera/$f.$v.derep.fasta;
    export fileOutUCDeNovo=$folderFilterChimera/$f.$v.uchime_denovo.uc;
    export fileOutChimerasDeNovo=$folderFilterChimera/$f.$v.uchime_denovo.chimeras.fasta;
    export fileOutNonchimerasDeNovo=$folderFilterChimera/$f.$v.uchime_denovo.nonchimeric.fasta;
    export fileOutNonchimeras=$folderSample/$f.$v.nonchimeric.fasta;
    
    #Run uchime_ref
    $USEARCH -uchime_ref $fileRAW -db $ChimeraReferenceUDB -threads 8 -strand plus -uchimeout $fileOutUCRef -chimeras $fileOutChimerasRef -nonchimeras $fileOutNonchimerasRef;
    #Run uchime_denovo
    #Dereplicate
    $USEARCH -derep_fulllength $fileRAW -uc $fileDerepOut -output $fileDerepFASTA -sizeout;
    #Filter chimeric sequences
    $USEARCH -uchime_denovo $fileDerepFASTA -uchimeout $fileOutUCDeNovo -chimeras $fileOutChimerasDeNovo -nonchimeras $fileOutNonchimerasDeNovo;
    
    #Get overlap between methods
    perl -e '
      #Read dereplication output file
      my %Rep;
      open (DEREP, $ENV{fileDerepOut}) or die;
      while (<DEREP>) {
        chomp; @fields = split /\t/;
        if ($fields[0] eq "H") {$Rep{$fields[9]}{$fields[8]} = 1}
      }
      close DEREP;
      
      #Read list of chimeras flagged by uchime_ref
      open (UCHIMEREF, $ENV{fileOutUCRef}) or die;
      while (<UCHIMEREF>) {chomp; @fields = split /\t/; if ($fields[17] eq "Y") {$Chimeras_ref{$fields[1]} = 1}}
      close UCHIMEREF;
      
      #Read list of chimeras flagged by uchime_denovo
      open (UCHIMEDENOVO, $ENV{fileOutUCDeNovo}) or die;
      while (<UCHIMEDENOVO>) {
        chomp; @fields = split /\t/;
        if ($fields[17] eq "Y") {
          $acc = $fields[1]; $acc =~ s/;.+//g;
          $Chimeras_denovo{$acc} = 1;
          foreach $rep (keys %{$Rep{$acc}}) {$Chimeras_denovo{$rep} = 1}
        }
      }
      close UCHIMEDENOVO;
      
      #Filter unaligned file using intersect of predicted chimeras
      open (UNALIGNED, $ENV{fileRAW}) or die;
      open (OUTPUT, ">$ENV{fileOutNonchimeras}") or die;
      while ($line = <UNALIGNED>) {
        chomp $line;
        if ($line =~ /^>/) {
          $line =~ s/^>//;
          #Look up whether chimera in both sets
          if (($Chimeras_denovo{$line} == 1) & ($Chimeras_ref{$line} == 1)) {next}
          $seqline = <UNALIGNED>;
          $length = length($seqline);
          #Skip if sequence is too short
          if ($length <= $ENV{minLength}) {next}
          print OUTPUT ">$line\n$seqline";
        }
      }
      close UNALIGNED;
      close OUTPUT;
    ';
    
    #Tidy up
    gzip $fileRAW;
    gzip $fileDerepOut;
    gzip $fileDerepFASTA;
    gzip $fileOutUCRef;
    gzip $fileOutUCDeNovo;
    gzip $fileOutChimerasRef;
    gzip $fileOutChimerasDeNovo;
    gzip $fileOutNonchimerasRef;
    gzip $fileOutNonchimerasDeNovo;
    gzip $fileOutNonchimeras;
  fi
done


