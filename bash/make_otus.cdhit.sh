#!/bin/bash
################################################################################
#Cluster sequences into OTUs using cd-hit, per sample
#=> at different thresholds
#
#Export (R-readable output files):
#=> sequence-to-OTU mapping
#=> OTU-to-sequence mapping
#
#
#2015-08-07
#sebastian.schmidt@imls.uzh.ch
################################################################################

f=$1;

#Set path to CD-HIT binary
CDHIT=/local/atlas/sebastian/bin/cd-hit/cd-hit-est;

#Set basic parameters
maxmem=100000;                    #maximum memory available (and yes, gaia has A LOT of memory)
wordsize=10;                      #kmer size for distance calculations
threads=8                         #number of threads for multicore cd-hit
lowestThresh=900                  #lowest threshold at which to output OTUs (x/1000)
highestThresh=1000                #highest threshold at which to output OTUs (x/1000)
increment=2                       #increment between thresholds at which OTUs are to be made (x/1000)

#Set basic directory and file names
folderBase=/mnt/mnemo3/sebastian/projects/tax.mapping/hmp;
folderData=$folderBase/1-data/per_sample;
folderSample=$folderData/$f;
folderOutput=$folderSample/cdhit;
mkdir $folderOutput;

#Iterate over subregions
for v in v13 v35 v69; do
  fileFASTAgz=$folderSample/$f.$v.nonchimeric.filtered.fasta.gz;
  if [ -e $fileSTOgz ]; then
    export V=$v;
    
    #Set file names
    export fileFASTA=$folderSample/$f.$v.tmp_cdhit.fasta;
    export SeqAccHash=$folderSample/$f.$v.nonchimeric.filtered.seq_acc.hash;
    export SeqIDHash=$folderSample/$f.$v.nonchimeric.filtered.seq_id.hash;
    
    #Gunzip
    gunzip -c $fileFASTAgz > $fileFASTA;
    
    #Iterate through cutoffs
    for c in $(seq $lowestThresh $increment $highestThresh); do
      #Format currnent cutoff
      cutoff=$(echo "$c/1000" | bc -l);
      
      #Set file names
      export fileOut=$folderOutput/$f.$v.cdhit.$c.fasta;
      export fileOutOTU=$folderOutput/$f.$v.cdhit.$c.otu;
      export fileOutMapping=$folderOutput/$f.$v.cdhit.$c.r.otu;
      export fileOutSeqMapping=$folderOutput/$f.$v.cdhit.$c.r.seq_mapping.otu;
      
      #Run CD-HIT
      $CDHIT -i $fileFASTA -o $fileOut -c $cutoff -T 8 -g 0 -M $maxmem -n $wordsize -d 0;
      
      #Process OTU data
      perl -e '
        #Load sequence data
        use Storable;
        %Acc = %{retrieve($ENV{SeqAccHash})};
        %ID = %{retrieve($ENV{SeqIDHash})};
        
        #Read CD-HIT "clstr" file
        open (CDHIT, $ENV{fileOut} . ".clstr") or die;
        open (OTU, ">$ENV{fileOutOTU}") or die;
        open (MAP, ">$ENV{fileOutMapping}") or die;
        $otu_idx = 0; $first = 1;
        while (<CDHIT>) {
          chomp;
          if (/^>/) {
            $otu_idx++;
            if ($first == 1) {print OTU ">OTU$otu_idx\n"; print MAP "$otu_idx\t"; $first = 0}
            else {print OTU "\n>OTU$otu_idx\n"; print MAP "\n$otu_idx\t"}
          }
          else {
            s/.+>//g; s/\.{3}.+//g;
            $acc = $_;
            $id = $Acc{$acc};
            print OTU "$acc\n";
            print MAP "$id,";
            $Seq{$id} = $otu_idx;
          }
        }
        close CDHIT;
        
        #Export sequence mapping
        open(SEQMAP, ">$ENV{fileOutSeqMapping}");
        foreach $id (sort {$a <=> $b} keys %Seq) {print SEQMAP "$id\t$Seq{$id}\n"}
        close SEQMAP;
      ';
      
      #Tidy up
      gzip -f $fileOut;
      gzip -f ${fileOut}.clstr;
      gzip -f ${fileOut}.bak.clstr;
      gzip -f $fileOutOTU;
    done
    
    #Tidy up
    rm $fileFASTA
  fi
done






