#!/bin/bash
################################################################################
#Global reference-based OTU demarcation using UCLUST
#=> at 97% (level "R")
#=> mapping against an UC- or AL-preclustered database of representative sequences
#
#Export:
#=> OTU files
#=> sequence-to-OTU mapping
#=> OTU-to-sequence mapping
#=> OTU datas
#=> OTU tables
#
#
#2015-11-24
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

#Turn "U" characters in databaset to "T", because UCLUST doesn't otherwise understand what it's supposed to do
#cat $fileDBtmp | perl -pi -e 's/U/T/g' > $fileDB;

#Iterate over subregions
for v in v13 v35; do
  export V=$v;
  
  #Set file names
  export fileSampleMapping=$folderData/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileUnaligned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.fasta;
  export fileUnalignedPruned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.pruned.fasta;
  export fileDerepFASTA=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.derep.fasta;
  export fileDerepCluster=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.derep.cluster;
  export fileTaxonomy=$folderOutput/$v.nonchimeric.filtered.full_db.taxonomy;
  export fileSeqAccHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_id.hash;
  
  #For V13 sequences, chop off leading ~70bp to obtain a sequence that UCLUST can match to the database
  if [ "$v" = "v13" ]; then
    cat $fileUnaligned | perl -ane 'if (/^>/) {print} else {chomp; $seq=$_; $l=length $seq; $start=$l-400; if($start<=0){print "$seq\n"} else{$keep=substr $seq,$start; print "$keep\n"}}' > $fileUnalignedPruned;
    #Dereplicate FASTA file
    $DEREPLICATE $fileUnalignedPruned $fileDerepFASTA $fileDerepCluster;
  elif [ "$v" = "v35" ]; then
    #Dereplicate FASTA file
    $DEREPLICATE $fileUnaligned $fileDerepFASTA $fileDerepCluster;
  fi
  
  #Split input file into parts of 20,000 sequences
  split --lines=40000 $fileDerepFASTA $fileDerepFASTA.split_
  
  for ref in uc_rep al_rep; do
    export REF=$ref;
    
    #Select current database
    if [ "$ref" = "uc_rep" ]; then
      export fileDB=$folderRefDB/full-ssu-db3.clustering.bacteria.uc.R.rep.fasta;
    elif [ "$ref" = "al_rep" ]; then
      export fileDB=$folderRefDB/full-ssu-db3.clustering.bacteria.al.R.rep.unaligned.U_to_T.fasta;
    fi
    
    #Set current file names
    export fileClosedUC=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.uc;
    export fileUnmappedFASTA=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.unmapped.fasta;
    export fileClosedRefOTU=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.otu;
    export fileClosedRefMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.r.otu;
    export fileClosedRefSeqMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.r.seq_mapping.otu;
    export fileClosedRefOTUTable=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.r.otu_table.tsv;
    export fileClosedRefOTUData=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.closed_ref.r.otu_data.tsv;
    
    #Run UCLUST, in parallel
    parallel -j20 "$UCLUST --input {} --id 0.97 --lib $fileDB --libonly --rev --maxaccepts $maxaccepts --maxrejects $maxrejects --w $wordlength --uc {}.$ref.uc" ::: $fileDerepFASTA.split_??
    #$UCLUST --input $fileDerepFASTA --id 0.97 --lib $fileDB --libonly --rev --maxaccepts $maxaccepts --maxrejects $maxrejects --w $wordlength --uc $fileClosedUC;
    
    #Merge *.uc output
    cat $fileDerepFASTA.split_??.$ref.uc > $fileClosedUC;
    
    #Process uc output
    #=> export OTU set, mappings, OTU table and OTU data
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
      @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
      open (TAX, $ENV{fileTaxonomy}) or die;
      while (<TAX>) {
        chomp;
        ($id, $acc, @fields) = split /\t/;
        foreach $tax_lev (@tax_levels) {$Seq{$acc}{Taxonomy}{$tax_lev} = shift @fields;}
      }
      close TAX;
      ##############################
      
      ##############################
      #Read dereplication clustering mapping
      open (CLUSTER, "$ENV{fileDerepCluster}") or die;
      while (<CLUSTER>) {
        chomp;
        ($seed, @children) = split /\t/;
        $Seed{$seed}{$seed} = 1;
        foreach $child (@children) {$Seed{$seed}{$child} = 1}
      }
      close CLUSTER;
      ##############################
      
      ##############################
      #Read dereplicated FASTA file
      open (DEREP, "$ENV{fileDerepFASTA}") or die;
      LINE:
      while ($line = <DEREP>) {
        chomp $line;
        if ($line =~ /^>/) {
          $line =~ s/^>//;
          $seqline = <DEREP>;
          chomp $seqline;
          $Seed_Seq{$line} = $seqline;
        }
      }
      close DEREP;
      ##############################
      
      ##############################
      #Read UC file & export unmapped FASTA file
      open (UC, $ENV{fileClosedUC}) or die;
      open (UNMAPPED, ">$ENV{fileUnmappedFASTA}") or die;
      LINE:
      while (<UC>) {
        chomp;
        next LINE if /^#/;
        @fields = split /\t/;
        if ($fields[0] eq "L") {push @otu_order => $fields[8] unless exists $seen_seeds{$fields[8]}; $seen_seeds{$fields[8]} = 1}
        elsif ($fields[0] eq "H") {foreach $child (keys %{$Seed{$fields[8]}}) {$raw_OTU{$fields[9]}{$child} = 1}}
        elsif ($fields[0] eq "N") {foreach $child (keys %{$Seed{$fields[8]}}) {print UNMAPPED ">$child\n$Seed_Seq{$fields[8]}\n"}}
      }
      close UC;
      close UNMAPPED;
      ##############################
      
      ##############################
      #Export OTU file, OTU mapping and collect sample mapping & taxonomy
      open(OTU, ">$ENV{fileClosedRefOTU}") or die;
      open(MAP, ">$ENV{fileClosedRefMapping}") or die;
      $otu_idx=1;
      foreach $seed (@otu_order) {
        $otu_name = $seed;
        #$prefix = "OTU_" . uc($ENV{V}) . "_UC_R_";
        #$otu_name = sprintf "$prefix%.6d", $otu_idx;
        
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
      open(SEQMAP, ">$ENV{fileClosedRefSeqMapping}") or die;
      foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
      close SEQMAP;
      ##############################
      
      ##############################
      #Export OTU data
      open(OTUDATA, ">$ENV{fileClosedRefOTUData}") or die;
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
      open(OTUTABLE, ">$ENV{fileClosedRefOTUTable}") or die;
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
    
    #Set file names for open reference
    export fileOpenUC=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.uc;
    export fileOpenRefOTU=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.otu;
    export fileOpenRefMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.r.otu;
    export fileOpenRefSeqMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.r.seq_mapping.otu;
    export fileOpenRefOTUTable=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.r.otu_table.tsv;
    export fileOpenRefOTUData=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.open_ref.r.otu_data.tsv;
    export fileUnmappedFASTAsorted=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_ref.R.$ref.unmapped.sorted.fasta;
    
    #Re-sort unmapped sequences by abundance
    $SORTABUND $fileUnmappedFASTA > $fileUnmappedFASTAsorted
    
    #Cluster unmapped sequences de novo
    $UCLUST --usersort --id 0.97 --input $fileUnmappedFASTAsorted --uc $fileOpenUC;
    
    #Process uc output
    #=> export OTU set, mappings, OTU table and OTU data
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
      @tax_levels = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species");
      open (TAX, $ENV{fileTaxonomy}) or die;
      while (<TAX>) {
        chomp;
        ($id, $acc, @fields) = split /\t/;
        foreach $tax_lev (@tax_levels) {$Seq{$acc}{Taxonomy}{$tax_lev} = shift @fields;}
      }
      close TAX;
      ##############################
      
      ##############################
      #Read dereplication clustering mapping
      open (CLUSTER, "$ENV{fileDerepCluster}") or die;
      while (<CLUSTER>) {
        chomp;
        ($seed, @children) = split /\t/;
        $Seed{$seed}{$seed} = 1;
        foreach $child (@children) {$Seed{$seed}{$child} = 1}
      }
      close CLUSTER;
      ##############################
      
      ##############################
      #Read dereplicated FASTA file
      open (DEREP, "$ENV{fileDerepFASTA}") or die;
      LINE:
      while ($line = <DEREP>) {
        chomp $line;
        if ($line =~ /^>/) {
          $line =~ s/^>//;
          $seqline = <DEREP>;
          chomp $seqline;
          $Seed_Seq{$line} = $seqline;
        }
      }
      close DEREP;
      ##############################
      
      ##############################
      #Read UC file for closed ref
      open (UC, $ENV{fileClosedUC}) or die;
      LINE:
      while (<UC>) {
        chomp;
        next LINE if /^#/;
        @fields = split /\t/;
        if ($fields[0] eq "L") {push @otu_order => $fields[8] unless exists $seen_seeds{$fields[8]}; $seen_seeds{$fields[8]} = 1}
        elsif ($fields[0] eq "H") {foreach $child (keys %{$Seed{$fields[8]}}) {$raw_OTU{$fields[9]}{$child} = 1}}
      }
      close UC;
      ##############################
      
      ##############################
      #Read UC file for open ref
      open (UC, $ENV{fileOpenUC}) or die;
      LINE:
      while (<UC>) {
        chomp;
        next LINE if /^#/;
        @fields = split /\t/;
        if ($fields[0] eq "S") {$raw_OTU{$fields[8]}{$fields[8]} = 1; push @otu_order => $fields[8]}
        elsif ($fields[0] eq "H") {$raw_OTU{$fields[9]}{$fields[8]} = 1}
      }
      close UC;
      ##############################
      
      ##############################
      #Export OTU file, OTU mapping and collect sample mapping & taxonomy
      open(OTU, ">$ENV{fileOpenRefOTU}") or die;
      open(MAP, ">$ENV{fileOpenRefMapping}") or die;
      $otu_idx=1;
      foreach $seed (@otu_order) {
        $otu_name = $seed;
        #$prefix = "OTU_" . uc($ENV{V}) . "_UC_R_";
        #$otu_name = sprintf "$prefix%.6d", $otu_idx;
        
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
      open(SEQMAP, ">$ENV{fileOpenRefSeqMapping}") or die;
      foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
      close SEQMAP;
      ##############################
      
      ##############################
      #Export OTU data
      open(OTUDATA, ">$ENV{fileOpenRefOTUData}") or die;
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
      open(OTUTABLE, ">$ENV{fileOpenRefOTUTable}") or die;
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
done









