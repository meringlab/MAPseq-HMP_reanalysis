#!/bin/bash
################################################################################
#Cluster sequences into OTUs using hpc-clust, globally
#=> at 97% similarity
#=> SL pre-clustering to 97%, then AL clustering of sub-partitions
#
#Export (R-readable output files):
#=> sequence-to-OTU mapping
#=> OTU-to-sequence mapping
#
#
#2015-08-18
#sebastian.schmidt@imls.uzh.ch
################################################################################

#Set paths to binaries
HPCCLUST=/local/erisdb/jfmrod/work/hpc-clust/hpc-clust;
FastTreeMP=/local/atlas/sebastian/bin/FastTreeMP;
makeNewick=/local/erisdb/jfmrod/work/microbial-atlas/make-hca-tree.sh;

#Set basic parameters
threshold=0.90                    #minimum threshold until which distances should be calculated
dfunc=gap                         #distance calculation function
threads=80                        #number of threads for multicore hpc-clust

#Set basic directory and file names
folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
folderData=$folderBase/1-data/per_sample;
folderOutput=$folderBase/1-data/global;

#Iterate over subregions
for v in v13 v35 v69; do
  #Set current file names
  export fileAlignment=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.fasta;
  export fileSampleMapping=$folderOutput/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileTaxonomy=$folderOutput/$v.nonchimeric.filtered.full_db.taxonomy;
  export fileSeqAccHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_id.hash;
  export fileOutClust=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust;
  export fileTMP=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.otu.tmp;
  export fileOTU=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.otu;
  export fileMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu;
  export fileSeqMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.seq_mapping.otu;
  export fileOTUTable=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu_table.tsv;
  export fileOTUData=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu_data.tsv;
  export fileRepTMP=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.rep.fasta.tmp;
  export fileRep=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.rep.fasta;
  export fileTree=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.tree;
  
  #Run hpc-clust
  $HPCCLUST -al true -cl true -sl true -t $threshold -dfunc $dfunc -nthreads $threads --ignoreMemThres 1 -ofile $fileOutClust $fileAlignment;
  
  #Make OTUs
  $HPCCLUST -makeotus $fileAlignment $fileOutClust.al 0.97 | grep -v "^#" | perl -pi -e 's/^>/\n>/g' | sed "1d" > $fileTMP;
  
  #Re-format current OTU file (=> OTU names)
  cat $fileTMP | ( perl -ane 'if (/^>/) {/>OTU(\d+)/; $idx=$1+1; printf ">OTU%.6d", $idx; print "\n"} else {print}'; echo ) > $fileOTU;
  rm $fileTMP;
  
  #Export OTU tables & co
  perl -e '
    use Storable;
    %Acc = %{retrieve($ENV{fileSeqAccHash})};
    %ID = %{retrieve($ENV{fileSeqIDHash})};
    
    #Get a fast subroutine for key lookup by max value
    sub largest_value (\%) {
        my $hash   = shift;
        my ($key, @keys) = keys   %$hash;
        my ($big, @vals) = values %$hash;
    
        for (0 .. $#keys) {
            if ($vals[$_] > $big) {
                $big = $vals[$_];
                $key = $keys[$_];
            }
        }
        $key
    }
    
    #Read sequence to sample mapping
    open(SAMPLEMAP, $ENV{fileSampleMapping}) or die;
    while (<SAMPLEMAP>) {chomp; ($acc, $sample) = split /\t/; $Seq{$acc}{Sample} = $sample}
    close SAMPLEMAP;
    
    #Read taxonomy per sequence
    @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
    open (TAX, $ENV{fileTaxonomy}) or die;
    while (<TAX>) {
      chomp;
      ($id, $acc, @fields) = split /\t/;
      foreach $tax_lev (@tax_levels) {$Seq{$acc}{Taxonomy}{$tax_lev} = shift @fields;}
    }
    close TAX;
    
    #Export mapping
    open(OTU, $ENV{fileOTU}) or die;
    open(MAP, ">$ENV{fileMapping}") or die;
    $otu_idx=1;
    LINE:
    while (<OTU>) {
      chomp;
      next LINE if /^#/;
      if (/^>/) {
        /^>(.+)/; $otu_name=$1;
        print MAP "$otu_idx\t";
        NEXTLINE:
        while ($acc = <OTU>) {
           chomp $acc;
           if ($acc =~ /^$/) {$OTU{$otu_name}{OTU_ID} = $otu_idx; $otu_idx++; print MAP "\n"; next LINE}
           $id = $Acc{$acc};
           print MAP "$id,";
           $Seq_to_OTU{$id} = $otu_idx;
           $sample = $Seq{$acc}{Sample};
           $Sample{$sample}{OTU}{$otu_name}++;
           $OTU{$otu_name}{Sample}{$sample}++;
           $OTU{$otu_name}{Size}++;
           foreach $tax_lev (@tax_levels) {$OTU{$otu_name}{Taxonomy}{$tax_lev}{$Seq{$acc}{Taxonomy}{$tax_lev}}++ unless ($Seq{$acc}{Taxonomy}{$tax_lev} ~~ ["", " "])}
        }
      }
    }
    close OTU; close MAP;
    
    #Export sequence mapping
    open(SEQMAP, ">$ENV{fileSeqMapping}") or die;
    foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
    close SEQMAP;
    
    #Export OTU data
    open(OTUDATA, ">$ENV{fileOTUData}") or die;
    print OTUDATA "OTU\tsize\tconsensus_domain\tconsensus_phylum\tconsensus_class\tconsensus_order\tconsensus_family\tconsensus_genus\tconsensus_species\tconsensus_taxonomy\n";
    foreach $otu_name (sort {$a cmp $b} keys %OTU) {
      print OTUDATA "$otu_name\t$OTU{$otu_name}{Size}";
      
      #Get consensus taxonomy
      $tax_string = "";
      foreach $tax_lev (@tax_levels) {
        if (exists $OTU{$otu_name}{Taxonomy}{$tax_lev}) {
          $largest = largest_value %{$OTU{$otu_name}{Taxonomy}{$tax_lev}};
          $dominance = $OTU{$otu_name}{Taxonomy}{$tax_lev}{$largest} / $OTU{$otu_name}{Size};
          
          if ($dominance > 0.5) {print OTUDATA "\t$largest"; $tax_string .= "$largest;"}
          else {print OTUDATA "\t"}
        } else {print OTUDATA "\t"}
      }
      $tax_string =~ s/;$//;
      print OTUDATA "\t$tax_string\n";
    }
    close OTUDATA;
    
    #Export OTU table
    open(OTUTABLE, ">$ENV{fileOTUTable}") or die;
    #Print header
    foreach $sample (sort {$a <=> $b} keys %Sample) {print OTUTABLE "\t$sample"; push @sample_order => $sample} print OTUTABLE "\n";
    foreach $otu (sort {$a cmp $b} keys %OTU) {
      print OTUTABLE "$otu";
      foreach $sample (@sample_order) {
        if (exists $OTU{$otu}{Sample}{$sample}) {print OTUTABLE "\t$OTU{$otu}{Sample}{$sample}"}
        else {print OTUTABLE "\t0"}
      }
      print OTUTABLE "\n";
    }
  ';
  
  #Get representative sequences
  $HPCCLUST -makereps $fileAlignment $fileOTU -nthreads $threads > $fileRepTMP;
  #Reformat FASTA file to harmonize OTU names
  perl -e '
    open (TMP, $ENV{fileRepTMP}) or die;
    $otu_idx = 1;
    while ($line = <TMP>) {
      chomp $line;
      next if $line =~ /^#/;
      next if $line =~ /^$/;
      if ($line =~ /^>/) {@fields=split / / => $line; shift @fields; printf ">OTU%.6d ", $otu_idx; print "\n"; $otu_idx++}
      else {print "$line\n"}
    }
  ' > $fileRep;
  
  #Get ML tree of representatives  
  $FastTreeMP -nt -gtr $fileRep > $fileTree;
  
  #Tidy up
  gzip -f $fileMapping;
  gzip -f $fileSeqMapping;
  gzip -f $fileOTUTable;
  gzip -f $fileOTU;
done





