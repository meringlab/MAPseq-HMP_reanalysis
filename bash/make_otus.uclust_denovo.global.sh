#!/bin/bash
################################################################################
#Global de novo clustering using UCLUST
#=> at 97% (level "R")
#=> matching the QIIME defaults as closely as possible
#
#Export:
#=> OTU files
#=> sequence-to-OTU mapping
#=> OTU-to-sequence mapping
#=> OTU datas
#=> OTU tables
#
#
#2015-08-27
#sebastian.schmidt@imls.uzh.ch
################################################################################

#Set path to uclust binary
UCLUST=~jfmrod/usr/bin/uclust;
SORTABUND=/local/erisdb/sebastian/toolbox_16S/2-perl/sort_fasta_by_abundance.pl;

#Set basic parameters to values used in QIIME HMP SOP
threshold=0.97                    #clustering threshold

#Set basic directory and file names
export folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
export folderData=$folderBase/1-data/global;
export folderOutput=$folderBase/1-data/global;

#Iterate over subregions
for v in v13 v35; do
  export V=$v;
  #Set current file names
  export fileSampleMapping=$folderData/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileUnaligned=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.fasta;
  export fileUnalignedSorted=$folderData/$v.nonchimeric.filtered.unaligned.global.denoised.sorted.fasta;
  export fileSampleMapping=$folderOutput/$v.nonchimeric.filtered.denoised.seq2sample.map;
  export fileTaxonomy=$folderOutput/$v.nonchimeric.filtered.full_db.taxonomy;
  export fileSeqAccHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_acc.hash;
  export fileSeqIDHash=$folderOutput/$v.nonchimeric.filtered.denoised.seq_id.hash;
  export fileUC=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.uc;
  export fileOTU=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.otu;
  export fileMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu;
  export fileSeqMapping=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.seq_mapping.otu;
  export fileOTUTable=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu_table.tsv;
  export fileOTUData=$folderOutput/$v.nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu_data.tsv;
  
  #Sort by abundance
  $SORTABUND $fileUnaligned > $fileUnalignedSorted;
  
  #Run UCLUST on sorted sequences
  $UCLUST --usersort --id 0.97 --input $fileUnalignedSorted --uc $fileUC;
  
  #Export OTU set, mappings, OTU table and OTU data
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
    #Read UC file
    open(UC, $ENV{fileUC}) or die;
    LINE:
    while (<UC>) {
      chomp;
      next LINE if /^#/;
      @fields = split /\t/;
      if ($fields[0] eq "S") {$raw_OTU{$fields[8]}{$fields[8]} = 1; push @otu_order => $fields[8]}
      else {$raw_OTU{$fields[9]}{$fields[8]} = 1}
    }
    close UC;
    ##############################
    
    ##############################
    #Export OTU file, OTU mapping and collect sample mapping & taxonomy
    open(OTU, ">$ENV{fileOTU}") or die;
    open(MAP, ">$ENV{fileMapping}") or die;
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
    open(SEQMAP, ">$ENV{fileSeqMapping}") or die;
    foreach $id (sort {$a <=> $b} keys %Seq_to_OTU) {print SEQMAP "$id\t$Seq_to_OTU{$id}\n"}
    close SEQMAP;
    ##############################
    
    ##############################
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
    ##############################
    
    ##############################
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
    ##############################
  ';
  
  
  ##############################
done




#Iterate over subregions
for v in v13 v35 v69; do
  fileFASTAgz=$folderSample/$f.$v.nonchimeric.filtered.fasta.gz;
  if [ -e $fileSTOgz ]; then
    export V=$v;
    
    #Set file names
    export fileFASTA=$folderSample/$f.$v.tmp_uclust_ref.fasta;
    export SeqAccHash=$folderSample/$f.$v.nonchimeric.filtered.seq_acc.hash;
    export SeqIDHash=$folderSample/$f.$v.nonchimeric.filtered.seq_id.hash;
    
    #Gunzip
    gunzip -c $fileFASTAgz > $fileFASTA;
    
    #Prune v13 sequences to work properly for uclust
    #=> cut away the ~30 first basepairs, to obtain a total length of 430bp
    if [ "$v" = "v13" ]; then
      cat $fileFASTA | perl -ane 'if (/^>/) {print} else {chomp; $seq=$_; $l=length $seq; $start=$l-430; if($start<=0){print "$seq\n"} else{$keep=substr $seq,$start; print "$keep\n"}}' > $folderSample/$f.$v.tmp_tmp.fasta;
      mv $folderSample/$f.$v.tmp_tmp.fasta $fileFASTA;
    fi
    
    #Reverse complement reads so that uclust can actually map them (I know, I know...)
    cat $fileFASTA | perl -pi -e 'unless (/^>/) {tr/ACGT/TGCA/; reverse}' > $folderSample/$f.$v.tmp_tmp.fasta;
    mv $folderSample/$f.$v.tmp_tmp.fasta $fileFASTA;
    
    #Iterate through reference databases
    for ref in uc al; do
      #Toggle reference database
      if [ "$ref" = "uc" ];
        then fileDB=$fileDBuc
      elif [ "$ref" = "al" ];
        then fileDB=$fileDBal;
      fi
      
      #Set current file names
      #Closed reference
      export fileClosedUC=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.uc;
      export fileClosedOTU=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.otu;
      export fileClosedOTUnames=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.names.otu;
      export fileClosedMapping=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.r.otu;
      export fileClosedMappingNames=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.r.names.otu;
      export fileClosedSeqMapping=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.seq_mapping.r.otu;
      export fileClosedSeqMappingNames=$folderOutput/$f.$v.uclust_ref.closed_ref.$mode.db_$ref.R.r.seq_mapping.names.otu;
      #Open reference
      export fileOpenUC=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.uc;
      export fileOpenOTU=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.otu;
      export fileOpenOTUnames=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.names.otu;
      export fileOpenMapping=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.r.otu;
      export fileOpenMappingNames=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.r.names.otu;
      export fileOpenSeqMapping=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.seq_mapping.r.otu;
      export fileOpenSeqMappingNames=$folderOutput/$f.$v.uclust_ref.open_ref.$mode.db_$ref.R.r.seq_mapping.names.otu;
      
      #Run UCLUST, closed reference
      $UCLUST --input $fileFASTA --id 0.97 --lib $fileDB --libonly --rev --maxaccepts $maxaccepts --maxrejects $maxrejects --w $wordlength --uc $fileClosedUC;
      
      #Process uc output
      perl -e '
        #Load sequence data
        use Storable;
        %Acc = %{retrieve($ENV{SeqAccHash})};
        %ID = %{retrieve($ENV{SeqIDHash})};
        
        #Read uc file
        open (UC, $ENV{fileClosedUC});
        while (<UC>) {
          chomp; @fields = split /\t/;
          if ($fields[0] eq "H") {$OTUs{$fields[9]}{$fields[8]} = 1}
        }
        close UC;
        
        #Export OTU files
        open (OTU, ">$ENV{fileClosedOTU}");
        open (OTUNAM, ">$ENV{fileClosedOTUnames}");
        open (MAP, ">$ENV{fileClosedMapping}");
        open (MAPNAM, ">$ENV{fileClosedMappingNames}");
        $otu_idx = 1;
        foreach $seed (keys %OTUs) {
          print OTU ">OTU$otu_idx\n";
          print OTUNAM ">$seed\n";
          print MAP "$otu_idx\t";
          print MAPNAM "$seed\t";
          foreach $acc (sort {$a cmp $b} keys %{$OTUs{$seed}}) {
            $id = $Acc{$acc};
            print OTU "$acc\n";
            print OTUNAM "$acc\n";
            print MAP "$id,";
            print MAPNAM "$id,";
            $Seq{$id}{OTU_Index} = $otu_idx;
            $Seq{$id}{OTU_Name} = $seed;
          }
          print OTU "\n";
          print OTUNAM "\n";
          print MAP "\n";
          print MAPNAM "\n";
          $otu_idx++;
        }
        close OTU;
        close OTUNAM;
        close MAP;
        close MAPNAM;
        
        #Export sequence mapping
        open(SEQMAP, ">$ENV{fileClosedSeqMapping}");
        open(SEQMAPNAM, ">$ENV{fileClosedSeqMappingNames}");
        foreach $id (sort {$a <=> $b} keys %Seq) {print SEQMAP "$id\t$Seq{$id}{OTU_Index}\n"; print SEQMAPNAM "$id\t$Seq{$id}{OTU_Name}\n"}
        close SEQMAP;
        close SEQMAPNAM;
      ';
      
      
      
      #Run UCLUST, open reference
      $UCLUST --input $fileFASTA --id 0.97 --rev --maxaccepts $maxaccepts --maxrejects $maxrejects --w $wordlength --uc $fileOpenUC;
      #Process uc output
      perl -e '
        #Load sequence data
        use Storable;
        %Acc = %{retrieve($ENV{SeqAccHash})};
        %ID = %{retrieve($ENV{SeqIDHash})};
        
        #Read uc file
        open (UC, $ENV{fileOpenUC});
        while (<UC>) {
          chomp; @fields = split /\t/;
          if ($fields[0] eq "H") {$OTUs{$fields[9]}{$fields[8]} = 1}
        }
        close UC;
        
        #Export OTU files
        open (OTU, ">$ENV{fileOpenOTU}");
        open (OTUNAM, ">$ENV{fileOpenOTUnames}");
        open (MAP, ">$ENV{fileOpenMapping}");
        open (MAPNAM, ">$ENV{fileOpenMappingNames}");
        $otu_idx = 1;
        foreach $seed (keys %OTUs) {
          print OTU ">OTU$otu_idx\n";
          print OTUNAM ">$seed\n";
          print MAP "$otu_idx\t";
          print MAPNAM "$seed\t";
          foreach $acc (sort {$a cmp $b} keys %{$OTUs{$seed}}) {
            $id = $Acc{$acc};
            print OTU "$acc\n";
            print OTUNAM "$acc\n";
            print MAP "$id,";
            print MAPNAM "$id,";
            $Seq{$id}{OTU_Index} = $otu_idx;
            $Seq{$id}{OTU_Name} = $seed;
          }
          print OTU "\n";
          print OTUNAM "\n";
          print MAP "\n";
          print MAPNAM "\n";
          $otu_idx++;
        }
        close OTU;
        close OTUNAM;
        close MAP;
        close MAPNAM;
        
        #Export sequence mapping
        open(SEQMAP, ">$ENV{fileOpenSeqMapping}");
        open(SEQMAPNAM, ">$ENV{fileOpenSeqMappingNames}");
        foreach $id (sort {$a <=> $b} keys %Seq) {print SEQMAP "$id\t$Seq{$id}{OTU_Index}\n"; print SEQMAPNAM "$id\t$Seq{$id}{OTU_Name}\n"}
        close SEQMAP;
        close SEQMAPNAM;
      ';
      
    done
    
    #Tidy up
    rm $fileFASTA;
  fi
done

