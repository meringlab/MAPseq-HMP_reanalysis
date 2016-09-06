#!/bin/bash
################################################################################
#Map sequences to pre-clustered reference using mapseq
#=> for multiple taxonomies/clusterings at the same time
#=> at different OTU thresholds at the same time
#=> map against full-ssu-db3, but also against uclust-preclustered database (97%)
#=> map to REPRESENTATIVE SEQUENCES ONLY (unlike complementary script on full databases)
#
#Export:
#=> mapping files at different thresholds
#=> taxonomic classification at different levels ("OTU data")
#=> FASTA files of unmapped sequences at 97% to be aligned in next step
#=> global, closed-reference OTU tables
#
#
#2015-11-19
#sebastian.schmidt@imls.uzh.ch
################################################################################


################################################################################
#Set path to binaries
MAPSEQ=~jfmrod/gaia/usr/bin/mapseq;
DENOISER=/home/jfmrod/gaia/usr/bin/denoiser;
HPCCLUST=/local/erisdb/jfmrod/work/hpc-clust/hpc-clust;
UCLUST=~jfmrod/usr/bin/uclust;
SORTABUND=/local/erisdb/sebastian/toolbox_16S/2-perl/sort_fasta_by_abundance.pl;

#Set parameters
export confThresh=0.6;
export consensusTaxThresh=0.5;
export maxMismatches=5;                       #eq to ca 1 allowed mismatch per 100bp
export minSeqFreq=3;
export dfunc=gap;

#Set basic file names
#Basic data folders
export folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
export folderData=$folderBase/1-data/global;
#Reference database folders & files
export folderRefDB=/local/erisdb/sebastian/toolbox_16S/reference_databases/full-ssu-db-ref;


#Make reference taxonomies of representative sequences per clustering method
#=> uc should have been done before
for cm in al cl sl; do
  export fileRep=$folderRefDB/full-ssu-db3.clustering.bacteria.$cm.R.rep.unaligned.fasta;
  export fileOTUTAX=$folderRefDB/full-ssu-db3.rep.$cm.R.otutax;
  echo "#cutoff: 0.97" > $fileOTUTAX;
  grep "^>" $fileRep | perl -ane 'chomp; s/^>//; print "$_\t$_\n"' >> $fileOTUTAX;
  
done

#Copy uc files from previous names, just to be consistent for subsequent looping
cp $folderRefDB/full-ssu-db3.otu_tax.uclust.R.rep.otutax $folderRefDB/full-ssu-db3.rep.uc.R.otutax;
cp $folderRefDB/full-ssu-db3.R.uclust.rep.fasta $folderRefDB/full-ssu-db3.clustering.bacteria.uc.R.rep.fasta;
################################################################################


################################################################################
#Process data
for v in v13 v35; do
  export V=$v;
  #Set current file names
  export fileUnaligned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.fasta;
  export fileSeqAccHash=$folderData/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderData/$v.nonchimeric.filtered.denoised.seq_id.hash;
  export fileAlignment=$folderData/$v.nonchimeric.filtered.aligned.global.raw.fasta;
  export fileUnalignedRAW=$folderData/$v.nonchimeric.filtered.unaligned.global.fasta;
  export fileSampleMapping=$folderData/$v.nonchimeric.filtered.denoised.seq2sample.map;
  
  #Iterate through clustering methods
  for cm in al cl sl uc; do
    export CM=$cm;
    #Set current file names
    export fileDB=$folderRefDB/full-ssu-db3.clustering.bacteria.$cm.R.rep.unaligned.fasta;
    export fileDBclust=$folderRefDB/full-ssu-db3.clustering.bacteria.$cm.R.rep.mscluster;
    export fileOTUTAX=$folderRefDB/full-ssu-db3.rep.$cm.R.otutax;
    export fileMSout=$folderData/$v.nonchimeric.filtered.rep.$cm.R.ms;
    
    #Run mapseq
    $MAPSEQ --nthreads 40 -cfile $fileDBclust $fileUnaligned $fileDB $fileOTUTAX > $fileMSout;
    
    
  done
  
  
  
  
  
  
  
  #Run mapseq
  $MAPSEQ --nthreads 40 -cfile $fileDBclust $fileUnaligned $fileDB $fileOTUTAX $fileOTUtaxAL $fileOTUtaxCL $fileOTUtaxSL $fileOTUtaxUC > $fileMSout;
  
  #Process output
  #=> export confident assignments per sequence, one file per taxonomy
  #=> export closed-ref OTU tables, one file per OTU taxonomy and threshold
  #=> export closed-ref OTU files, one file per OTU taxonomy and threshold
  #=> export closed-ref OTU data files, one file per OTU taxonomy and threshold
  #=> export closed-ref sequence-to-OTU and OTU-to-sequence mappings, one file per OTU taxonomy and threshold
  #=> export (non-denoised) alignments for unmapped sequences, one file per OTU taxonomy and threshold
  perl -e '
    #Get a fast subroutine for key lookup by max value (for lazy people like me)
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
    
    #Load sequence accession and ID hashes
    use Storable;
    %Acc = %{retrieve($ENV{fileSeqAccHash})};
    %ID = %{retrieve($ENV{fileSeqIDHash})};
    
    ##############################
    #Pre-define taxonomic levels
    @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
    @otu_levels = ("F", "G", "R", "S", "T");
    ##############################
    
    ##############################
    #Load full sequence alignment
    open (FASTA, $ENV{fileAlignment}) or die;
    while ($acc=<FASTA>) {chomp $acc; if ($acc =~ /^>/) {$acc =~ s/^>//; $seq=<FASTA>; chomp $seq; $Seq{$acc}{Sequence}=$seq}}
    close FASTA;
    ##############################
    
    ##############################
    #Load full unaligned sequences file
    open (FASTA, $ENV{fileUnalignedRAW}) or die;
    while ($acc=<FASTA>) {chomp $acc; if ($acc =~ /^>/) {$acc =~ s/^>//; $seq=<FASTA>; chomp $seq; $Seq{$acc}{Unaligned}=$seq}}
    close FASTA;
    ##############################
    
    ##############################
    #Load full sequence to sample mapping
    open (SAMPLEMAP, $ENV{fileSampleMapping}) or die;
    while (<SAMPLEMAP>) {chomp; ($acc, $sample) = split /\t/; $Seq{$acc}{Sample} = $sample}
    close SAMPLEMAP;
    ##############################
    
    ##############################
    #Set current file names
    $file_seq_tax_tax = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.taxonomy";
    $file_seq_tax_al = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.al.tax";
    $file_seq_tax_cl = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.cl.tax";
    $file_seq_tax_sl = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.sl.tax";
    $file_seq_tax_uc = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.uc.tax";
    ##############################
    #Load and process mapseq output
    ##############################
    #Open files
    open (MAPSEQ, $ENV{fileMSout}) or die;
    open (SEQTAXTAX, ">$file_seq_tax_tax") or die;
    open (SEQTAXAL, ">$file_seq_tax_al") or die;
    open (SEQTAXCL, ">$file_seq_tax_cl") or die;
    open (SEQTAXSL, ">$file_seq_tax_sl") or die;
    open (SEQTAXUC, ">$file_seq_tax_uc") or die;
    ##############################
    #Open files to export unmapped sequences
    for $cm ("al", "cl", "sl") {
      for $otu_lev (@otu_levels) {
        $curr_file = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.$cm.$otu_lev.full_db.unmapped.fasta";
        open ($files_unmapped{$cm}{$otu_lev}, ">$curr_file") or die;
      }
    }
    $curr_file = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.uc.R.full_db.unmapped.fasta";
    open ($files_unmapped{uc}{R}, ">$curr_file") or die;
    ##############################
    #Loop through all sequences
    while (<MAPSEQ>) {
      chomp;
      @tax_fields = split /\t\t/;
      
      #Extract current query sequence accession
      $acc = (split /\t/ => $tax_fields[0])[0]; $id = $Acc{$acc};
      foreach $fh (SEQTAXTAX, SEQTAXAL, SEQTAXCL, SEQTAXSL, SEQTAXUC) {print $fh "$id\t$acc"}
      
      #Get sample for current sequence
      $sample = $Seq{$acc}{Sample};
      
      ##############################
      #Process taxonomy
      @fields = split /\t/ => $tax_fields[1];
      %seq_tax = ();
      #Iterate over taxonomic levels and export
      $export_string = "";
      foreach $tax_lev (@tax_levels) {
        #Get assigned taxonomy
        $tl = shift @fields;
        #Get confidence of assignment
        $conf = shift @fields;
        #Throw away fields which are not used 
        $pred_thresh = shift @fields; $dash = shift @fields;
        
        #Check if assignment was above confidence threshold
        if ($conf > $ENV{confThresh}) {$export_string .= "$tl;"; print SEQTAXTAX "\t$tl"; $Seq{$acc}{Taxonomy}{$tax_lev} = $tl; $seq_tax{$tax_lev} = $tl}
        else {$export_string .= ";"; $seq_tax{$tax_lev} = "NA"; print SEQTAXTAX "\t"}
      }
      $export_string =~ s/;+$//;
      $Seq{$acc}{Taxonomy}{Consensus} = $export_string;
      print SEQTAXTAX "\t$export_string\n";
      ##############################
      #Process AL, CL and SL OTU mapping
      CM:
      foreach $cm ("al", "cl", "sl") {
        #Get current OTU taxonomy mapping results fields
        if ($cm eq "al") {@fields = split /\t/ => $tax_fields[2]; $fh = SEQTAXAL}
        elsif ($cm eq "cl") {@fields = split /\t/ => $tax_fields[3]; $fh = SEQTAXCL}
        elsif ($cm eq "sl") {@fields = split /\t/ => $tax_fields[4]; $fh = SEQTAXSL}
        #Iterate over taxonomic levels and export
        $export_string = "";
        foreach $otu_lev (@otu_levels) {
          #Get assigned taxonomy
          $otu = shift @fields;
          #Get confidence of assignment
          $conf = shift @fields;
          #Throw away fields which are not used 
          $pred_thresh = shift @fields; $dash = shift @fields;
          
          #Check if assignment was above confidence threshold
          if ($conf > $ENV{confThresh}) {
            $export_string .= "$otu;";
            $OTU{$cm}{$otu_lev}{$otu}{Sequences}{$acc} = 1;
            $OTU{$cm}{$otu_lev}{$otu}{Samples}{$sample}++;
            $Sample{$sample}{OTU}{$cm}{$otu_lev}{$otu}++;
            for $tax_lev (@tax_levels) {$OTU{$cm}{$otu_lev}{$otu}{Taxonomy}{$tax_lev}{$seq_tax{$tax_lev}}++ unless $seq_tax{$tax_lev} eq "NA"}
            print $fh "\t$otu";
          } else {
            $export_string .= ";";
            print $fh "\t";
            $fh_unmapped = $files_unmapped{$cm}{$otu_lev};
            print $fh_unmapped ">$acc\n$Seq{$acc}{Sequence}\n";
          }
        }
        $export_string =~ s/;+$//;
        print $fh "\t$export_string\n";
      }
      ##############################
      #Process UC OTU mapping
      @fields = split /\t/ => $tax_fields[5];
      #Get assigned taxonomy
      $otu = shift @fields;
      #Get confidence of assignment
      $conf = shift @fields;
      if ($conf > $ENV{confThresh}) {
        $OTU{uc}{R}{$otu}{Sequences}{$acc} = 1;
        $OTU{uc}{R}{$otu}{Samples}{$sample}++;
        $Sample{$sample}{OTU}{uc}{R}{$otu}++;
        for $tax_lev (@tax_levels) {$OTU{uc}{R}{$otu}{Taxonomy}{$tax_lev}{$seq_tax{$tax_lev}}++ unless $seq_tax{$tax_lev} eq "NA"}
        print SEQTAXUC "\t$otu\n";
      } else {
        print SEQTAXUC "\t\n";
        $fh_unmapped = $files_unmapped{uc}{R};
        print $fh_unmapped ">$acc\n$Seq{$acc}{Unaligned}\n";
      }
    }
    ##############################
    close MAPSEQ;
    close SEQTAXTAX; close SEQTAXAL; close SEQTAXCL; close SEQTAXSL; close SEQTAXUC;
    foreach $cm (keys %files_unmapped) {foreach $otu_lev (keys %{$files_unmapped{$cm}}) {close $files_unmapped{$cm}{$otu_lev}}}
    ##############################
    
    
    ##############################
    #Export OTU data for AL, CL, SL & UC
    CM:
    foreach $cm ("al", "cl", "sl", "uc") {
      #Iterate through OTU thresholds
      for $otu_lev (sort {$a cmp $b} keys %{$OTU{$cm}}) {
        #Set file names
        $file_otu = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.$cm.$otu_lev.full_db.closed_ref.otu";
        $file_otu_mapping = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.$cm.$otu_lev.full_db.closed_ref.mapping.otu";
        $file_otu_data = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.$cm.$otu_lev.full_db.closed_ref.otu_data.tsv";
        $file_otu_table = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.$cm.$otu_lev.full_db.closed_ref.otu_table.tsv";
        #Open files
        open (OTU, ">$file_otu") or die;
        open (OTUMAP, ">$file_otu_mapping") or die;
        open (OTUDATA, ">$file_otu_data") or die;
        open (OTUTABLE, ">$file_otu_table") or die;
        
        #Print OTU data header
        print OTUDATA "OTU\tsize\tconsensus_domain\tconsensus_phylum\tconsensus_class\tconsensus_order\tconsensus_family\tconsensus_genus\tconsensus_species\tconsensus_taxonomy\n";
        
        #Print OTU table header
        foreach $sample (sort {$a <=> $b} keys %Sample) {print OTUTABLE "\t$sample"; push @sample_order => $sample} print OTUTABLE "\n";
        
        #Iterate through matched OTUs
        foreach $otu (sort {$a cmp $b} keys %{$OTU{$cm}{$otu_lev}}) {
          print OTU ">$otu\n";
          print OTUMAP "$otu\t";
          print OTUDATA "$otu";
          print OTUTABLE "$otu";
          
          #Export...
          #=> to OTU file
          #=> to OTU mapping file
          $otu_size = 0;
          foreach $acc (sort {$a cmp $b} keys %{$OTU{$cm}{$otu_lev}{$otu}{Sequences}}) {
            print OTU "$acc\n";
            print OTUMAP "$Acc{$acc},";
            $otu_size++;
          }
          print OTU "\n";
          print OTUMAP "\n";
          print OTUDATA "\t$otu_size";
          
          #Get consensus taxonomy
          $tax_string = "";
          foreach $tax_lev (@tax_levels) {
            if (exists $OTU{$cm}{$otu_lev}{$otu}{Taxonomy}{$tax_lev}) {
              $largest = largest_value %{$OTU{$cm}{$otu_lev}{$otu}{Taxonomy}{$tax_lev}};
              $dominance = $OTU{$cm}{$otu_lev}{$otu}{Taxonomy}{$tax_lev}{$largest} / $otu_size;
              
              if ($dominance > $ENV{consensusTaxThresh}) {print OTUDATA "\t$largest"; $tax_string .= "$largest;"}
              else {print OTUDATA "\t"}
            } else {print OTUDATA "\t"}
          }
          $tax_string =~ s/;$//;
          print OTUDATA "\t$tax_string\n";
          
          #Export to OTU table
          foreach $sample (@sample_order) {
            if (exists $OTU{$cm}{$otu_lev}{$otu}{Samples}{$sample}) {print OTUTABLE "\t$OTU{$cm}{$otu_lev}{$otu}{Samples}{$sample}"}
            else {print OTUTABLE "0\t"}
          }
          print OTUTABLE "\n";
        }
        close OTU; close OTUMAP; close OTUDATA; close OTUTABLE;
      }
    }
    
    ##############################
  ';
  
  
  #Generate open reference OTU sets
  #=> for R level only (that is the one needed for benchmarking)
  #=> per AL, CL & SL
  for cm in al cl sl; do
    export CM=$cm;
    export fileUnmappedRAW=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.unmapped.fasta;
    export fileDenoised=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.denoised.cluster;
    export fileAlignmentDenoised=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.unmapped.aligned.denoised.fasta;
    export fileOutClust=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.unmapped.aligned.denoised;
    export fileTMP=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.unmapped.aligned.denoised.$cm.tmp;
    export fileClosedRefOTU=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.closed_ref.otu;
    export fileOTU=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.open_ref.otu;
    export fileOTUmapping=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.open_ref.mapping.otu;
    export fileOTUseqmapping=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.open_ref.seq_mapping.otu;
    export fileOTUtable=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.open_ref.otu_table.tsv;
    export fileOTUdata=$folderData/$v.nonchimeric.filtered.$cm.R.full_db.open_ref.otu_data.tsv;
    
    #Run denoiser to de-replicate alignment
    $DENOISER -aligned true -maxmiss $maxMismatches -nthreads 40 -of tmp $fileUnmappedRAW > $fileDenoised;
    
    #Process de-replication output to generate updated alignment
    perl -e '
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
      open(FASTA, $ENV{fileUnmappedRAW}) or die;
      while ($acc = <FASTA>) {
        chomp $acc;
        if ($acc =~ /^>/) {$acc =~ s/^>//; if (exists $Seed{$acc}) {$seq=<FASTA>; chomp $seq; $Seed{$acc}{Seq}=$seq}}
      }
      close FASTA;
      
      #Export sequences & list of retained sequences
      open(OUT, ">$ENV{fileAlignmentDenoised}") or die;
      foreach $acc (@seq_order) {print OUT ">$acc\n$Seed{$Child{$acc}}{Seq}\n"}
      close OUT;
    ';
    
    #Run HPC-Clust on processed alignments
    $HPCCLUST -$cm true -t 0.97 -dfunc $dfunc -nthreads 40 --ignoreMemThres 1 -ofile $fileOutClust $fileAlignmentDenoised;
    
    #Make OTUs
    $HPCCLUST -makeotus $fileAlignmentDenoised $fileOutClust.al 0.97 | grep -v "^#" | perl -pi -e 's/^>/\n>/g' | sed "1d" > $fileTMP;
    
    #Re-format current OTU file (=> OTU names)
    cat $fileClosedRefOTU > $fileOTU;
    cat $fileTMP | ( perl -ane 'if (/^>/) {/>OTU(\d+)/; $idx=$1+1; $prefix="OTU_" . uc($ENV{V}) . "_" . uc($ENV{CM}) . "_R_"; printf ">$prefix%.6d", $idx; print "\n"} else {print}'; echo ) >> $fileOTU;
    rm $fileTMP;
    
    #Export open reference OTU set, mappings, OTU table and OTU data
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
      
      ##############################
      #Read sequence to sample mapping
      open(SAMPLEMAP, $ENV{fileSampleMapping}) or die;
      while (<SAMPLEMAP>) {chomp; ($acc, $sample) = split /\t/; $Seq{$acc}{Sample} = $sample}
      close SAMPLEMAP;
      ##############################
      
      ##############################
      #Read taxonomy per sequence
      $file_seq_tax = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.taxonomy";
      @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
      open (TAX, $file_seq_tax) or die;
      while (<TAX>) {
        chomp;
        ($id, $acc, @fields) = split /\t/;
        foreach $tax_lev (@tax_levels) {$Seq{$acc}{Taxonomy}{$tax_lev} = shift @fields;}
      }
      close TAX;
      ##############################
      
      ##############################
      #Export OTU mapping and collect sample mapping & taxonomy
      open(OTU, $ENV{fileOTU}) or die;
      open(MAP, ">$ENV{fileOTUmapping}") or die;
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
             next NEXTLINE if $id ~~ ["", " "];
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
      ##############################
      
      ##############################
      #Export sequence mapping
      open(SEQMAP, ">$ENV{fileOTUseqmapping}") or die;
      foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
      close SEQMAP;
      ##############################
      
      ##############################
      #Export OTU data
      open(OTUDATA, ">$ENV{fileOTUdata}") or die;
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
      ##############################
      
      ##############################
      #Export OTU table
      open(OTUTABLE, ">$ENV{fileOTUtable}") or die;
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
      ##############################
    ';
  done
  
  #Generate open reference OTU set for UCLUST
  #=> R level / 97% only
  export fileUnmappedRAW=$folderData/$v.nonchimeric.filtered.uc.R.full_db.unmapped.fasta;
  export fileUnmappedSorted=$folderData/$v.nonchimeric.filtered.uc.R.full_db.unmapped.sorted.fasta;  
  export fileUC=$folderData/$v.nonchimeric.filtered.uc.R.full_db.unmapped.uc;
  export fileClosedRefOTU=$folderData/$v.nonchimeric.filtered.uc.R.full_db.closed_ref.otu;
  export fileOTU=$folderData/$v.nonchimeric.filtered.uc.R.full_db.open_ref.otu;
  export fileOTUmapping=$folderData/$v.nonchimeric.filtered.uc.R.full_db.open_ref.mapping.otu;
  export fileOTUseqmapping=$folderData/$v.nonchimeric.filtered.uc.R.full_db.open_ref.seq_mapping.otu;
  export fileOTUtable=$folderData/$v.nonchimeric.filtered.uc.R.full_db.open_ref.otu_table.tsv;
  export fileOTUdata=$folderData/$v.nonchimeric.filtered.uc.R.full_db.open_ref.otu_data.tsv;
  
  #Sort sequences by abundance
  $SORTABUND $fileUnmappedRAW > $fileUnmappedSorted;
  
  #Run UCLUST on remaining sequences
  $UCLUST --usersort --id 0.97 --input $fileUnmappedSorted --uc $fileUC;
  
  #Export open reference OTU set, mappings, OTU table and OTU data
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
    
    ##############################
    #Read sequence to sample mapping
    open(SAMPLEMAP, $ENV{fileSampleMapping}) or die;
    while (<SAMPLEMAP>) {chomp; ($acc, $sample) = split /\t/; $Seq{$acc}{Sample} = $sample}
    close SAMPLEMAP;
    ##############################
    
    ##############################
    #Read taxonomy per sequence
    $file_seq_tax = "$ENV{folderData}/$ENV{V}.nonchimeric.filtered.full_db.taxonomy";
    @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
    open (TAX, $file_seq_tax) or die;
    while (<TAX>) {
      chomp;
      ($id, $acc, @fields) = split /\t/;
      foreach $tax_lev (@tax_levels) {$Seq{$acc}{Taxonomy}{$tax_lev} = shift @fields;}
    }
    close TAX;
    ##############################
    
    ##############################
    #Read UC file
    open(UC, $ENV{fileUC}) or die;
    LINE:
    while (<UC>) {
      chomp;
      next LINE if /^#/;
      @fields = split /\t/;
      if ($fields[0] eq "S") {$raw_OTU{$fields[8]}{$fields[8]} = 1; push @otu_order => $fields[8]}
      else {$raw_OTU{$fields[8]}{$fields[9]} = 1}
    }
    ##############################
    
    ##############################
    #Export OTU file, OTU mapping and collect sample mapping & taxonomy
    open(OTU, ">$ENV{fileOTU}") or die;
    open(MAP, ">$ENV{fileOTUmapping}") or die;
    $otu_idx=1;
    foreach $seed (@otu_order) {
      $prefix = "OTU_" . uc($ENV{V}) . "_UC_R_";
      $otu_name = sprintf "$prefix%.6d", $otu_idx;
      
      print OTU ">$otu_name\n";
      print MAP "$otu_idx\t";
      ACC:
      foreach $acc (sort {$a cmp $b} keys %{$raw_OTU{$seed}}) {
        $id = $Acc{$acc};
        next ACC if $id ~~ ["", " "];
        print OTU "$acc\n";
        print MAP "$id,";
        
        $OTU{$otu_name}{OTU_ID} = $otu_idx;
        $Seq_to_OTU{$id} = $otu_idx;
        $sample = $Seq{$acc}{Sample};
        $Sample{$sample}{OTU}{$otu_name}++;
        $OTU{$otu_name}{Sample}{$sample}++;
        $OTU{$otu_name}{Size}++;
        foreach $tax_lev (@tax_levels) {$OTU{$otu_name}{Taxonomy}{$tax_lev}{$Seq{$acc}{Taxonomy}{$tax_lev}}++ unless ($Seq{$acc}{Taxonomy}{$tax_lev} ~~ ["", " "])}
      }
      print OTU "\n";
      print MAP "\n";
      
      $otu_idx++;
    }
    close OTU; close MAP;
    ##############################
    
    ##############################
    #Export sequence mapping
    open(SEQMAP, ">$ENV{fileOTUseqmapping}") or die;
    foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
    close SEQMAP;
    ##############################
    
    ##############################
    #Export OTU data
    open(OTUDATA, ">$ENV{fileOTUdata}") or die;
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
    ##############################
    
    ##############################
    #Export OTU table
    open(OTUTABLE, ">$ENV{fileOTUtable}") or die;
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
    ##############################
  ';
done








################################################################################









