#!/usr/bin/perl
################################################################################
#Process HMP raw data
#=> start out w/ "HM16STR" data -> pre-processed FASTA files, resolved per samples
#=> sort per-sample FASTA files by sequenced subregion
#=> remove (potentially) mislabelled or misclassified samples
#
#Export:
#=> sample metadata table
#=> data in one folder per sample, subfolders per subregion
#=> lists of samples per subregion
#=> global mapping file of sequences per sample
#
#2015-08-04
#sebastian.schmidt@imls.uzh.ch
################################################################################


################################################################################
#Begin
use strict;
use warnings;
no warnings 'uninitialized';
use List::Util qw(min max sum shuffle);
use List::MoreUtils;
use POSIX qw(strftime);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError) ;
$| = 1;
################################################################################


################################################################################
#Set basic paths & file names
my $folderBase = "/mnt/mnemo3/sebastian/projects/tax.mapping/hmp";
my $folderRAW = $folderBase . "/0-raw";
my $folderParam = $folderBase . "/parameter_files";
my $folderData = $folderBase . "/1-data/per_sample";

#Set basic parameters
my $min_length = 200;                     #this corresponds to the HMP defaults (mothur & QIIME pipelines)

#Preallocate global data structures
my %Sample_Data = ();
my %Sequence_Data = ();
my %contaminated_samples = ();
my %list_SRS = ();
my %list_samples = ();
################################################################################


################################################################################
#Read DACC parameter files

#Get list of mislabelled or contaminated samples
my $file_contaminated = $folderParam . "/mislabeled_contaminated_list_combined.txt";
open (FILE, $file_contaminated) or die;
while (my $line = <FILE>) {chomp $line; next if $line =~ /^#/; my ($srs, undef) = split /\./ => $line; $contaminated_samples{$srs} = 1}
close FILE;

#Get a list of available samples (SRS accessions)
opendir(DIR, $folderRAW);
while (my $fn = readdir(DIR)) {chomp $fn; next unless $fn =~ /\.gz$/; $fn =~ s/\.fsa\.gz//; $list_SRS{$fn} = 1}
closedir(DIR);

#Read sample metadata
#v35 mapping file
my $file_ppAll_V35 = $folderParam . "/ppAll_V35_map.tsv";
open(FILE, $file_ppAll_V35) or die;
while (my $line = <FILE>) {
  chomp $line;
  next if $line =~ /^Sample/;
  my @fields = split /\t/ => $line;
  #Skip if there is no sequence data available for current SRS
  next unless ($list_SRS{$fields[6]} == 1);
  
  #Store data
  my $sn = $fields[3];
  $Sample_Data{$sn}{RSID} = $fields[1];
  $Sample_Data{$sn}{PSN} = $fields[2];
  $Sample_Data{$sn}{Experiment_Accession}{$fields[4]} = 1;
  $Sample_Data{$sn}{Run_ID}{$fields[5]} = 1;
  $Sample_Data{$sn}{Sex} = $fields[10];
  $Sample_Data{$sn}{Body_Subsite} = $fields[11];
  $Sample_Data{$sn}{Body_Site} = $fields[12];
  $Sample_Data{$sn}{Visit} = $fields[13];
}
close FILE;

#v13 mapping file => only updates (samples that are not in v35)
my $file_ppAll_V13 = $folderParam . "/ppAll_V13_map.tsv";
open(FILE, $file_ppAll_V13) or die;
while (my $line = <FILE>) {
  chomp $line;
  next if $line =~ /^Sample/;
  my @fields = split /\t/ => $line;
  #Skip if there is no sequence data available for current SRS
  next unless ($list_SRS{$fields[6]} == 1);
  #Skip if this sample was already annotated in v35
  next if exists ($Sample_Data{$fields[3]});
  
  #Store data
  my $sn = $fields[3];
  $Sample_Data{$sn}{RSID} = $fields[1];
  $Sample_Data{$sn}{PSN} = $fields[2];
  $Sample_Data{$sn}{Experiment_Accession}{$fields[4]} = 1;
  $Sample_Data{$sn}{Run_ID}{$fields[5]} = 1;
  $Sample_Data{$sn}{Sex} = $fields[10];
  $Sample_Data{$sn}{Body_Subsite} = $fields[11];
  $Sample_Data{$sn}{Body_Site} = $fields[12];
  $Sample_Data{$sn}{Visit} = $fields[13];
}
close FILE;
################################################################################


################################################################################
#Process per-sample FASTA files

#Iterate through all raw *.gz files
opendir(DIR, $folderRAW);
while (my $tmp_fn = readdir(DIR)) {
  chomp $tmp_fn; next unless $tmp_fn =~ /\.gz$/;
  my $fn = "$folderRAW/$tmp_fn";
  my $fh = $fn; $fh =~ s/\.gz$//;
  
  print "\n>Processing $fn...\n";
  
  #Gunzip
  gunzip $fn => $fh or die "gunzip failed: $GunzipError\n";
  #Open file
  open(FILE, $fh) or die;
  #Read data
  LINE:
  while (my $line = <FILE>) {
    if ($line =~ /^>/) {
      #Get sequence line
      my $seq = <FILE>; chomp $seq;
      #Skip if sequence is too short
      next LINE if length($seq) < $min_length;
      #Get current data
      my ($acc, undef, $sn, undef, $primer, undef) = split / / => $line;
      $sn =~ s/^sample=//; $primer =~ s/^primer=//;
      #Skip if current SN has no metadata
      unless (exists $Sample_Data{$sn}) {print "Unknown SN for $acc:\t$sn\n"; next LINE}
      #Extract V-region
      my $v;
      if ($primer eq "V1-V3") {$v = "v13"; $list_samples{$v}{$sn}++}
      elsif ($primer eq "V3-V5") {$v = "v35"; $list_samples{$v}{$sn}++}
      elsif ($primer eq "V6-V9") {$v = "v69"; $list_samples{$v}{$sn}++}
      else {print "Unknown primer set for $acc:\t$primer\n"; next LINE}
      #Store
      $Sample_Data{$sn}{Sequences}{$v}{$acc} = 1;
      $Sequence_Data{$acc}{SN} = $sn;
      $Sequence_Data{$acc}{V} = $v;
      $Sequence_Data{$acc}{Seq} = $seq;
    }
  }
  close FILE;
  
  #Gzip again
  gzip $fh => $fn or die "gzip failed: $GzipError\n";
  
  print "done.\n";
}
closedir(DIR);

#Export list of samples w/ v13, v35 and v69 data
my @v = ("v13", "v35", "v69");
foreach my $v (@v) {
  open(FILE, ">$folderParam/sample_list.$v") or die;
  foreach my $sn (sort {$a cmp $b} keys %{$list_samples{$v}}) {print FILE "$sn\t$list_samples{$v}{$sn}\n"}
  close FILE;
}

#Iterate through samples and export FASTA files per sample & subregion
SN:
foreach my $sn (sort {$a cmp $b} keys %Sample_Data) {
  #Skip if current sample has no sequences
  unless (exists $Sample_Data{$sn}{Sequences}) {print "Skipping $sn:\t no sequences found.\n"; next SN}
  print "\nExporting $sn\n";
  #Make target directory
  my $curr_dir = "$folderData/$sn";
  mkdir $curr_dir;
  #Iterate through subregions
  V:
  foreach my $v (sort {$a cmp $b} keys %{$Sample_Data{$sn}{Sequences}}) {
    print "$v:\t$list_samples{$v}{$sn}\t"; my $count = 0;
    my $curr_file = "$curr_dir/$v.fasta";
    open(FILE, ">$curr_file") or die;
    ACC:
    foreach my $acc (sort {$a cmp $b} keys %{$Sample_Data{$sn}{Sequences}{$v}}) {print FILE "$acc\n$Sequence_Data{$acc}{Seq}\n"; $count++}
    close FILE;
    print "$count\n";
  }
}

#Iterate through all sequence accessions (globally) and export sample mapping
print "\nExporting sequence-to-sample mapping...\n";
open(FILE, ">$folderParam/sample_mapping.raw.tsv") or die;
foreach my $acc (sort {$a cmp $b} keys %Sequence_Data) {my $tmp_acc = $acc; $tmp_acc =~ s/^>//; print FILE "$tmp_acc\t$Sequence_Data{$acc}{V}\t$Sequence_Data{$acc}{SN}\n"}
close FILE;
print "done.\n\n";

#Export sample metadata
print "\nExporting sample metadata...\n";
open(FILE, ">$folderParam/sample_data.raw.tsv") or die;
#Print header
print FILE "SN\t" .
  "PSN\t" .
  "RSID\t" .
  "Sex\t" .
  "Visit\t" .
  "Body_Site\t" .
  "Body_Subsite\t" .
  "Raw_v13\t" .
  "Raw_v35\t" .
  "Raw_v69\n";
SN:
foreach my $sn (sort {$a cmp $b} keys %Sample_Data) {
  print FILE "$sn\t" .
    "$Sample_Data{$sn}{PSN}\t" .
    "$Sample_Data{$sn}{RSID}\t" .
    "$Sample_Data{$sn}{Sex}\t" .
    "$Sample_Data{$sn}{Visit}\t" .
    "$Sample_Data{$sn}{Body_Site}\t" .
    "$Sample_Data{$sn}{Body_Subsite}\t";
  
  if (exists $list_samples{v13}{$sn}) {print FILE "$list_samples{v13}{$sn}\t"} else {print FILE "NA\t"}
  if (exists $list_samples{v35}{$sn}) {print FILE "$list_samples{v35}{$sn}\t"} else {print FILE "NA\t"}
  if (exists $list_samples{v69}{$sn}) {print FILE "$list_samples{v69}{$sn}\n"} else {print FILE "NA\n"}
}
close FILE;
print "done.\n\n";


################################################################################







