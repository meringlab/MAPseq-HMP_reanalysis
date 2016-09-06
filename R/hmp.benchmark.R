#!/usr/bin/Rscript
################################################################################
################################################################################
#Benchmark of different OTU demarcation methods on the HMP dataset
#
#The following methods are tested (all at a clustering threshold of 97%):
#=> hpc-clust, de novo, average linkage
#=> uclust, de novo
#=> mapseq, closed reference, AL-reference, full db
#=> mapseq, closed reference, UC-reference, full db
#=> mapseq, open reference, AL-reference, full db
#=> mapseq, open reference, UC-reference, full db (clustering of unmapped sequences using uclust)
#=> uclust, closed reference, AL-reference, representative sequences db
#=> uclust, closed reference, UC-reference, representative sequences db
#=> uclust, open reference, AL-reference, representative sequences db (clustering of unmapped sequences using uclust)
#=> uclust, open reference, UC-reference, representative sequences db
#
#Tests are done on HMP samples, comparing datasets on the V13 and V35 subregions of the 16S SSU rRNA
#The following analyses are performed:
#=> fraction of total (confidently) mapped sequences per sample for V13 and V35 (closed ref methods only)
#=> agreement in OTUs and OTU abundance per sample between V13 and V35 (closed ref methods only)
#=> agreement of alpha diversity estimates per sample between V13 and V35 (all methods)
#=> agreement of beta diversity estimates per sample pair between V13 and V35 (all methods)
#=> agreement of partitions (AMI) between de novo, closed ref and open ref per sample and V region (triplets of methods)
#
#
#2016-02-05
#sebastian.schmidt@imls.uzh.ch
################################################################################
################################################################################


################################################################################
################################################################################
# Load Packages
library("foreach");
library("doMC");
library("parallel");
library("Matrix");
library("ggplot2");
library("grid");
library("gridExtra");
library("RColorBrewer");
library("data.table");
library("vegan");
library("reshape2");
################################################################################
################################################################################


################################################################################
################################################################################
#Set parameters
PARAM <- list();
PARAM$folder.input <- "/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/1-data/global/";
PARAM$folder.results <- "/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/5-results/";
PARAM$file.sample_metadata <- "/mnt/mnemo3/sebastian/projects/tax.mapping/hmp/parameter_files/sample_data.raw.tsv";
PARAM$thresh.sample_size <- 1000;
PARAM$use.cores <- 40;

################################################################################
################################################################################


################################################################################
################################################################################
#Define functions

split.feature.vector <- function (X) {eval(parse(text=paste("c(", substr(X, 1, nchar(X)-1), ")")))}

#Chao 1 richness estimator
chao.1 <- function(tmp.freq) {
  #Remove OTUs that have not been observed (freq == 0)
  freq <- tmp.freq[tmp.freq > 0];
  s.obs <- length(freq);
  n.singl <- length(which(freq == 1));
  n.doubl <- length(which(freq == 2));
  s.chao <- s.obs + ((n.singl * (n.singl - 1)) / (2 * (n.doubl + 1)));
  return(s.chao)
}

#Simpson
simpson <- function(tmp.freq) {
  #Remove OTUs that have not been observed (freq == 0)
  freq <- tmp.freq[tmp.freq > 0];
  s.obs <- length(freq);
  n.tot <- sum(freq);
  s.simpson <- sum(freq * (freq - 1)) / (n.tot * (n.tot -1));
  s.simpson
}

#Shannon
shannon <- function(tmp.freq) {
  #Remove OTUs that have not been observed (freq == 0)
  freq <- tmp.freq[tmp.freq > 0];
  s.obs <- length(freq);
  n.tot <- sum(freq);
  n.rel.tmp <- freq / n.tot; n.rel <- n.rel.tmp[n.rel.tmp > 0];
  s.shannon <- -(sum(n.rel * log(n.rel)));
}

#Hill diversity ("true diversity")
hill.diversity <- function(tmp.freq, q) {
  if (q == 1) {warning("Hill diversity is not defined for q = 1. Calculating the limit case of D = exp(H), the exponential of the Shannon entropy, instead."); return(exp(shannon(tmp.freq)))}
  #Remove OTUs that have not been observed (freq == 0)
  freq <- tmp.freq[tmp.freq > 0];
  n.tot <- sum(freq);
  qD <- sum((freq / n.tot) ^ q) ^ (1 / (1-q));
  qD
}

################################################################################
################################################################################


################################################################################
################################################################################
#Load data
#=> sample metadata
#=> OTU tables
#=> OTU data
################################################################################
################################################################################
#Sample metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, sep="\t");
rownames(sample.data.raw) <- paste0("SN_", as.character(sample.data.raw$SN));

###################################
#OTU data
#=> per method
data <- list(); data[["v13"]] <- list(); data[["v35"]] <- list();
t <- 1;
###################################
#De novo, average linkage
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "de_novo.al";
  data[[sr]][[t]]$mode <- "de_novo";
  data[[sr]][[t]]$tool <- "hpc_clust";
  data[[sr]][[t]]$reference <- "none";
  data[[sr]][[t]]$ref_scope <- "none";
  data[[sr]][[t]]$algorithm <- "average_linkage";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu_table.tsv.gz"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.otu"), what=list(integer(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
  #Sequence mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.hpc_clust.0_97.r.seq_mapping.otu"), what=list(integer(), integer()), sep="\t");
  data[[sr]][[t]]$seq_mapping <- tmp.data[[2]];
  data[[sr]][[t]]$seq_mapping[data[[sr]][[t]]$seq_mapping == 0] <- NA;
}
t <- t + 1;
#Get total number of sequences per subregion
n.seq.v13 <- sum(colSums(data[["v13"]][[1]]$ot));
n.seq.v35 <- sum(colSums(data[["v35"]][[1]]$ot));
###################################
#De novo, uclust
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "de_novo.uc";
  data[[sr]][[t]]$mode <- "de_novo";
  data[[sr]][[t]]$tool <- "uclust";
  data[[sr]][[t]]$reference <- "none";
  data[[sr]][[t]]$ref_scope <- "none";
  data[[sr]][[t]]$algorithm <- "uclust";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.otu"), what=list(integer(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
  #Sequence mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_denovo.R.r.seq_mapping.otu"), what=list(integer(), integer()), sep="\t");
  data[[sr]][[t]]$seq_mapping <- tmp.data[[2]];
  data[[sr]][[t]]$seq_mapping[data[[sr]][[t]]$seq_mapping == 0] <- NA;
}
t <- t + 1;
###################################
#Closed reference, mapseq, al reference, full db
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "closed_ref.mapseq.al.full_db";
  data[[sr]][[t]]$mode <- "closed_ref";
  data[[sr]][[t]]$tool <- "mapseq";
  data[[sr]][[t]]$reference <- "al";
  data[[sr]][[t]]$ref_scope <- "full_db";
  data[[sr]][[t]]$algorithm <- "NA";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.closed_ref.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.closed_ref.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.closed_ref.mapping.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Closed reference, mapseq, uclust reference, full db
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "closed_ref.mapseq.uc.full_db";
  data[[sr]][[t]]$mode <- "closed_ref";
  data[[sr]][[t]]$tool <- "mapseq";
  data[[sr]][[t]]$reference <- "uclust";
  data[[sr]][[t]]$ref_scope <- "full_db";
  data[[sr]][[t]]$algorithm <- "NA";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.closed_ref.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.closed_ref.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.closed_ref.mapping.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Closed reference, uclust, al reference, representatives
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "closed_ref.uclust.al.rep";
  data[[sr]][[t]]$mode <- "closed_ref";
  data[[sr]][[t]]$tool <- "uclust";
  data[[sr]][[t]]$reference <- "al";
  data[[sr]][[t]]$ref_scope <- "rep";
  data[[sr]][[t]]$algorithm <- "NA";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.closed_ref.r.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.closed_ref.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.closed_ref.r.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Closed reference, uclust, uclust reference, representatives
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "closed_ref.uclust.uclust.rep";
  data[[sr]][[t]]$mode <- "closed_ref";
  data[[sr]][[t]]$tool <- "uclust";
  data[[sr]][[t]]$reference <- "uclust";
  data[[sr]][[t]]$ref_scope <- "rep";
  data[[sr]][[t]]$algorithm <- "NA";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.closed_ref.r.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.closed_ref.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.closed_ref.r.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Open reference, mapseq, al reference, full db
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "open_ref.mapseq.al.full_db";
  data[[sr]][[t]]$mode <- "open_ref";
  data[[sr]][[t]]$tool <- "mapseq";
  data[[sr]][[t]]$reference <- "al";
  data[[sr]][[t]]$ref_scope <- "full_db";
  data[[sr]][[t]]$algorithm <- "al";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.open_ref.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.open_ref.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.al.R.full_db.open_ref.mapping.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Open reference, mapseq, uclust reference, full db
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "open_ref.mapseq.uc.full_db";
  data[[sr]][[t]]$mode <- "open_ref";
  data[[sr]][[t]]$tool <- "mapseq";
  data[[sr]][[t]]$reference <- "uclust";
  data[[sr]][[t]]$ref_scope <- "full_db";
  data[[sr]][[t]]$algorithm <- "uclust";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.open_ref.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.open_ref.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.uc.R.full_db.open_ref.mapping.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Open reference, uclust, al reference, representatives
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "open_ref.uclust.al.rep";
  data[[sr]][[t]]$mode <- "open_ref";
  data[[sr]][[t]]$tool <- "uclust";
  data[[sr]][[t]]$reference <- "al";
  data[[sr]][[t]]$ref_scope <- "rep";
  data[[sr]][[t]]$algorithm <- "uclust";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.open_ref.r.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.open_ref.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.al_rep.open_ref.r.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
t <- t + 1;
###################################
#Open reference, uclust, uclust reference, representatives
for (sr in c("v13", "v35")) {
  data[[sr]][[t]] <- list();
  data[[sr]][[t]]$v <- sr;
  data[[sr]][[t]]$name <- "open_ref.uclust.uclust.rep";
  data[[sr]][[t]]$mode <- "open_ref";
  data[[sr]][[t]]$tool <- "uclust";
  data[[sr]][[t]]$reference <- "uclust";
  data[[sr]][[t]]$ref_scope <- "rep";
  data[[sr]][[t]]$algorithm <- "uclust";
  #OTU table
  data[[sr]][[t]]$ot <- Matrix(as.matrix(read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.open_ref.r.otu_table.tsv"), header=T, check.names=F, sep="\t", row.names=1)), sparse=T);
  colnames(data[[sr]][[t]]$ot) <- paste0("SN_", colnames(data[[sr]][[t]]$ot));
  #OTU data
  data[[sr]][[t]]$otu_data <- read.table(file = paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.open_ref.r.otu_data.tsv"), header=F, skip=1, sep="\t", row.names=1);
  #OTU mapping
  tmp.data <- scan(file=paste0(PARAM$folder.input, sr, ".nonchimeric.filtered.aligned.global.denoised.uc_ref.R.uc_rep.open_ref.r.otu"), what=list(character(), character()), sep="\t");
  data[[sr]][[t]]$otu_mapping <- mclapply(tmp.data[[2]], function(str) {as.numeric(unlist(strsplit(str, ",", fixed=T)))}, mc.cores = 40);
}
n.datasets <- t;
###################################

###################################
#Filter data
#=> identify samples represented in both v13 and v35 at sufficient depths (>1,000 sequences)
#=> prune sample metadata table accordingly
###################################
tmp.sums.v13 <- colSums(data[["v13"]][[1]]$ot);
tmp.sums.v35 <- colSums(data[["v35"]][[1]]$ot);
use.samples <- sort(intersect(names(tmp.sums.v13[tmp.sums.v13 > PARAM$thresh.sample_size]), names(tmp.sums.v35[tmp.sums.v35 > PARAM$thresh.sample_size])));
data.sample <- sample.data.raw[use.samples, ];
###################################
#Save
save(data, data.sample, PARAM, file=paste0(PARAM$folder.results, "all_data.RData"));
################################################################################
################################################################################


################################################################################
################################################################################
#Plot total fractions of (confidently) mapped sequences per sample
#=> for closed reference methods only
#=> relative to the full "filtered" sample size (as used in de novo and open ref clustering)
#=> per subregion
################################################################################
################################################################################
#Preallocate data.frame to collect results
res.frac_mapped <- as.data.frame(matrix(nrow=0, ncol=ncol(data.sample)+8));
colnames(res.frac_mapped) <- c(colnames(data.sample), "v", "name", "mode", "tool", "reference", "ref_scope", "algorithm", "fraction_mapped");
###################################
#Iterate through OTU datasets and collect fractions of mapped sequences
for (sr in c("v13", "v35")) {
  #Get reference sample sizes (number of sequences retained after all filtering steps)
  curr.ref.sizes <- colSums(data[[sr]][[1]]$ot[, use.samples]);
  for (i in seq(1, n.datasets)) {
    #Process only closed reference datasets
    if (data[[sr]][[i]]$mode != "closed_ref") {next()}
    #Calculate fraction of mapped reads per sample
    curr.frac <- colSums(data[[sr]][[i]]$ot[, use.samples]) / curr.ref.sizes;
    #Store
    res.frac_mapped <- rbind(res.frac_mapped, data.frame(
      data.sample,
      v=data[[sr]][[i]]$v,
      name=data[[sr]][[i]]$name,
      mode="closed_ref",
      tool=data[[sr]][[i]]$tool,
      reference=data[[sr]][[i]]$reference,
      ref_scope=data[[sr]][[i]]$ref_scope,
      algorithm=data[[sr]][[i]]$algorithm,
      fraction_mapped=curr.frac)
    );
  }
}
###################################
#Plot
curr.plot <- ggplot(res.frac_mapped, aes(x=v, y=fraction_mapped, fill=name, color=name)) +
  geom_boxplot(alpha=0.7, outlier.colour=NA) +
  geom_point(size=2, alpha=0.5, position=position_jitterdodge()) +
  xlab("V region") +
  ylab("Fraction of confidently mapped sequences") +
  theme_bw();
ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste0(PARAM$folder.results, "fraction_of_mapped_sequences.pdf"), useDingbats=FALSE);
###################################
#Reshape and export as *.tsv file
res.frac_mapped.wide <- dcast(res.frac_mapped, SN + PSN + RSID + Sex + Visit + Body_Site + Body_Subsite + Raw_v13 + Raw_v35 ~ interaction(v, name), value.var="fraction_mapped");
write.table(res.frac_mapped.wide, file=paste0(PARAM$folder.results, "fraction_of_mapped_sequences.tsv"), col.names=NA, sep = "\t", quote=F);
###################################
#Save
save(res.frac_mapped, file=paste0(PARAM$folder.results, "fraction_of_mapped_sequences.RData"));
################################################################################
################################################################################


################################################################################
################################################################################
#Test agreement between subregions per sample
#=> for closed reference methods only
#=> per sample, compare OTUs in v13 and v35 (mapped to the same common reference)
#=> get correlation of OTU relative abundances per sample
#=> get overlap as Jaccard & Bray-Curtis indices
################################################################################
################################################################################
#Preallocate data.frame to collect results
res.overlap <- as.data.frame(matrix(nrow=0, ncol=ncol(data.sample)+10));
colnames(res.overlap) <- c(colnames(data.sample), "v", "name", "mode", "tool", "reference", "ref_scope", "algorithm", "abundance_correlation", "jaccard", "bray_curtis");
sample.sizes.v13 <- colSums(data[["v13"]][[1]]$ot[, use.samples]);
sample.sizes.v35 <- colSums(data[["v35"]][[1]]$ot[, use.samples]);
###################################
#Iterate through OTU datasets
for (i in seq(1, n.datasets)) {
  #Process only closed reference datasets
  if (data[["v13"]][[i]]$mode != "closed_ref") {next()}
  #Extract OTU tables
  ot.v13 <- data[["v13"]][[i]]$ot[, use.samples];
  ot.v35 <- data[["v35"]][[i]]$ot[, use.samples];
  union.otus <- union(rownames(ot.v13), rownames(ot.v35));
  tmp.v13 <- tmp.v35 <- numeric(length(union.otus));
  names(tmp.v13) <- names(tmp.v35) <- union.otus;
  #Normalize counts by global (filtered) sample sizes
  ot_rel.v13 <- t(t(ot.v13) / sample.sizes.v13);
  ot_rel.v35 <- t(t(ot.v35) / sample.sizes.v35);
  #Correlate relative abundances per sample
  curr.cor <- unlist(lapply(use.samples, function(s) {
    d.v13 <- tmp.v13; d.v13[rownames(ot.v13)] <- ot_rel.v13[, s];
    d.v35 <- tmp.v35; d.v35[rownames(ot.v35)] <- ot_rel.v35[, s];
    cor(d.v13, d.v35, method="pearson")
  }));
  #Get Jaccard overlap per sample
  curr.jac <- unlist(lapply(use.samples, function(s) {
    length(intersect(rownames(ot.v13)[ot.v13[,s] > 0], rownames(ot.v35)[ot.v35[,s] > 0])) / length(union(rownames(ot.v13)[ot.v13[,s] > 0], rownames(ot.v35)[ot.v35[,s] > 0]))
  }));
  #Get Bray-Curtis similarity per sample
  curr.bc <- unlist(lapply(use.samples, function(s) {
    a.N <- sample.sizes.v13[s]; b.N <- sample.sizes.v35[s];
    a.freq <- tmp.v13; a.freq[rownames(ot.v13)] <- ot.v13[, s];
    b.freq <- tmp.v35; b.freq[rownames(ot.v35)] <- ot.v35[, s];
    2 * (sum(pmin(a.freq, b.freq)) / (a.N + b.N))
  }));
  #Store data
  res.overlap <- rbind(res.overlap, data.frame(
    data.sample,
    v="both",
    name=data[["v13"]][[i]]$name,
    mode="closed_ref",
    tool=data[["v13"]][[i]]$tool,
    reference=data[["v13"]][[i]]$reference,
    ref_scope=data[["v13"]][[i]]$ref_scope,
    algorithm=data[["v13"]][[i]]$algorithm,
    abundance_correlation=curr.cor,
    jaccard=curr.jac,
    bray_curtis=curr.bc
  ));
}
###################################
#Plot
###################################
#Abundance correlations
curr.plot <- ggplot(res.overlap, aes(x=tool, y=abundance_correlation, fill=name, color=name)) +
  geom_boxplot(alpha=0.7, outlier.colour=NA) +
  geom_point(size=2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0.3, position="dodge") +
  xlab("Mapping tool") +
  ylab("Abundance correlation of mapped OTUs between v13 and v35") +
  theme_bw();
ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste0(PARAM$folder.results, "overlap.abundance_correlation.pdf"), useDingbats=FALSE);
###################################
#Jaccard overlap
curr.plot <- ggplot(res.overlap, aes(x=tool, y=jaccard, fill=name, color=name)) +
  geom_boxplot(alpha=0.7, outlier.colour=NA) +
  geom_point(size=2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0.3, position="dodge") +
  xlab("Mapping tool") +
  ylab("Jaccard overlap of mapped OTUs between v13 and v35") +
  theme_bw();
ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste0(PARAM$folder.results, "overlap.jaccard.pdf"), useDingbats=FALSE);
###################################
#Bray-Curtis
curr.plot <- ggplot(res.overlap, aes(x=tool, y=bray_curtis, fill=name, color=name)) +
  geom_boxplot(alpha=0.7, outlier.colour=NA) +
  geom_point(size=2, alpha=0.5, position=position_jitterdodge()) +
  geom_violin(alpha=0.3, position="dodge") +
  xlab("Mapping tool") +
  ylab("Bray-Curtis similarity of mapped OTUs between v13 and v35") +
  theme_bw();
ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste0(PARAM$folder.results, "overlap.bray_curtis.pdf"), useDingbats=FALSE);
###################################
#Save
save(res.overlap, file=paste0(PARAM$folder.results, "overlap_between_subregions.RData"));
################################################################################
################################################################################


################################################################################
################################################################################
#Test agreement of alpha diversity estimates per sample between subregions
#=> for all methods
#=> calculate N_obs, Chao.1, Shannon and Simpson per sample
#=> correlate alpha div between v13 and v35 per sample
################################################################################
################################################################################
#Preallocate results and data collectors
data.alpha_div <- as.data.frame(matrix(nrow=0, ncol=ncol(data.sample)+15));
colnames(data.alpha_div) <- c(colnames(data.sample), "v", "name", "mode", "tool", "reference", "ref_scope", "algorithm", "v13.n_obs", "v35.n_obs", "v13.chao1", "v35.chao1", "v13.shannon", "v35.shannon", "v13.hill_2", "v35.hill_2");
res.cor.alpha_div <- as.data.frame(matrix(nrow=n.datasets, ncol=11));
colnames(res.cor.alpha_div) <- c("v", "name", "mode", "tool", "reference", "ref_scope", "algorithm", "n_obs", "chao1", "shannon", "hill_2");
###################################
#Iterate through OTU datasets
for (i in seq(1, n.datasets)) {
  #Preallocate current results collector
  curr.alpha_div <- data.frame(
    data.sample,
    v="both",
    name=data[["v13"]][[i]]$name,
    mode=data[["v13"]][[i]]$mode,
    tool=data[["v13"]][[i]]$tool,
    reference=data[["v13"]][[i]]$reference,
    ref_scope=data[["v13"]][[i]]$ref_scope,
    algorithm=data[["v13"]][[i]]$algorithm,
    v13.n_obs=NA, v35.n_obs=NA,
    v13.chao1=NA, v35.chao1=NA,
    v13.shannon=NA, v35.shannon=NA,
    v13.hill_2=NA, v35.hill_2=NA
  )
  #Iterate over subregions
  for (sr in c("v13", "v35")) {
    #Number of observed OTUs (N_obs)
    curr.alpha_div[use.samples, paste0(sr, ".n_obs")] <- apply(data[[sr]][[i]]$ot[, use.samples], 2, function(ov) {length(which(ov > 0))});
    #Chao1 richness estimator
    curr.alpha_div[use.samples, paste0(sr, ".chao1")] <- apply(data[[sr]][[i]]$ot[, use.samples], 2, chao.1);
    #Shannon entropy
    curr.sample.sizes <- eval(parse(text=paste0("sample.sizes.", sr)));
    curr.alpha_div[use.samples, paste0(sr, ".shannon")] <- sapply(use.samples, function(s) {
      tmp.freq <- data[[sr]][[i]]$ot[, s];
      freq <- tmp.freq[tmp.freq > 0];
      s.obs <- length(freq);
      n.tot <- curr.sample.sizes[s];
      n.rel.tmp <- freq / n.tot; n.rel <- n.rel.tmp[n.rel.tmp > 0];
      -(sum(n.rel * log(n.rel)));
    });
    #Hill diversity, order 2
    curr.alpha_div[use.samples, paste0(sr, ".hill_2")] <- sapply(use.samples, function(s) {
      tmp.freq <- data[[sr]][[i]]$ot[, s];
      freq <- tmp.freq[tmp.freq > 0];
      n.tot <- curr.sample.sizes[s];
      sum((freq / n.tot) ^ 2) ^ (1 / (1-2))
    });
  }
  #Store in global alpha diversity collector
  data.alpha_div <- rbind(data.alpha_div, curr.alpha_div);
  #Correlate diversity estimates per sample between subregions
  res.cor.alpha_div[i, 1:7] <- c("both", data[["v13"]][[i]]$name, data[["v13"]][[i]]$mode, data[["v13"]][[i]]$tool, data[["v13"]][[i]]$reference, data[["v13"]][[i]]$ref_scope, data[["v13"]][[i]]$algorithm);
  #N_obs
  res.cor.alpha_div[i, "n_obs"] <- cor(curr.alpha_div[, "v13.n_obs"], curr.alpha_div[, "v35.n_obs"], method="pearson");
  #Chao1
  res.cor.alpha_div[i, "chao1"] <- cor(curr.alpha_div[, "v13.chao1"], curr.alpha_div[, "v35.chao1"], method="pearson");
  #Shannon
  res.cor.alpha_div[i, "shannon"] <- cor(curr.alpha_div[, "v13.shannon"], curr.alpha_div[, "v35.shannon"], method="pearson");
  #Hill, q=2
  res.cor.alpha_div[i, "hill_2"] <- cor(curr.alpha_div[, "v13.hill_2"], curr.alpha_div[, "v35.hill_2"], method="pearson");
}
###################################
#Calculate log2(FC) of alpha diversity estimators between samples
#=> log2(v13 / v35)
data.alpha_div$log2FC.n_obs <- log2(data.alpha_div$v13.n_obs / data.alpha_div$v35.n_obs);
data.alpha_div$log2FC.chao1 <- log2(data.alpha_div$v13.chao1 / data.alpha_div$v35.chao1);
data.alpha_div$log2FC.shannon <- log2(data.alpha_div$v13.shannon / data.alpha_div$v35.shannon);
data.alpha_div$log2FC.hill_2 <- log2(data.alpha_div$v13.hill_2 / data.alpha_div$v35.hill_2);
#Plot
for (alpha in c("n_obs", "chao1", "shannon", "hill_2")) {
  curr.plot <- ggplot(data.alpha_div, aes_string(x="mode", y=paste0("log2FC.", alpha), fill="name", color="name")) +
    geom_boxplot(alpha=0.7, outlier.colour=NA) +
    geom_point(size=2, alpha=0.5, position=position_jitterdodge()) +
    xlab("OTU demarcation mode") +
    ylab(paste0("Log2(v13/v35) for ", alpha)) +
    theme_bw();
  ggsave(plot=curr.plot, height=10, width=20, dpi=300, filename=paste0(PARAM$folder.results, "alpha_div.log2FC.", alpha, ".pdf"), useDingbats=FALSE);
}
###################################
#Plot alpha-div scatters, in a large grid
curr.plots <- list();
#Collect plots
jj <- 1;
#Iterate through OTU datasets
for (i in seq(1, n.datasets)) {
  name=data[["v13"]][[i]]$name;
  #Iterate through alpha diversity indices
  for (alpha in c("n_obs", "chao1", "shannon", "hill_2")) {
    if (alpha == "n_obs") {limits <- c(0,5)}
    if (alpha == "chao1") {limits <- c(0,5)}
    if (alpha == "shannon") {limits <- c(-3, 1)}
    if (alpha == "hill_2") {limits <- c(0, 7)}
    curr.data <- data.frame(v13=log10(data.alpha_div[data.alpha_div$name == name, paste0("v13.", alpha)]), v35=log10(data.alpha_div[data.alpha_div$name == name, paste0("v35.", alpha)]));
    my_grob = grobTree(textGrob(as.character(res.cor.alpha_div[i, alpha]), x=0.1,  y=0.95, hjust=0, gp=gpar(fontsize=15)));
    curr.plots[[jj]] <- ggplot(curr.data, aes(x=v13, y=v35)) + geom_point(alpha=0.4) + scale_colour_hue(l=50) + scale_y_continuous(limits=limits) + scale_x_continuous(limits=limits) + geom_smooth(method=lm) + ggtitle(paste(name, alpha)) + annotation_custom(my_grob) + theme_bw();
    jj <- jj+1;
  }
}
#Export plot
pdf(file = paste0(PARAM$folder.results, "alpha_div.scatter.pdf"), width=60, height=100);
do.call(grid.arrange, c(curr.plots, nrow=n.datasets, ncol=4));
dev.off();
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate partition similarities
#=> Adjusted Mutual Information (AMI) of closed and open ref vs de novo
################################################################################
################################################################################
#Preallocate results collector
res.ami <- matrix(nrow=2, ncol=4);
rownames(res.ami) <- c("average_linkage", "uclust");
colnames(res.ami) <- c("closed_full", "open_full", "closed_rep", "open_rep");
res.ami.list <- list();
#Preallocate comparisons to perform
make.comparisons <- matrix(nrow=8, ncol=2);
make.comparisons[1,] <- c(1,3);
make.comparisons[2,] <- c(2,4);
make.comparisons[3,] <- c(1,7);
make.comparisons[4,] <- c(2,8);
make.comparisons[5,] <- c(1,5);
make.comparisons[6,] <- c(2,6);
make.comparisons[7,] <- c(1,9);
make.comparisons[8,] <- c(2,10);
#Iterate through subregions
for (sr in c("v13", "v35")) {
  res.ami.list[[sr]] <- list(); res.ami.list[[sr]]$res <- res.ami;
  #Iterate through comparisons
  for (i in seq(1, nrow(make.comparisons))) {
    #Get current OTU mapping (for set 2, ref-based clustering)
    curr.otu_map <- data[[sr]][[make.comparisons[i,2]]]$otu_mapping;
    #Get OTU sizes for set 2
    sizes.2 <- unlist(lapply(curr.otu_map, length));
    #Get current sequence mapping (for set 1, de novo clustering)
    curr.seq_map <- data[[sr]][[make.comparisons[i,1]]]$seq_mapping;
    #Get OTU sizes for set 1
    #=> corrected to take into account only sequences encountered in set 2
    sizes.1 <- tabulate(curr.seq_map[unlist(curr.otu_map)]);
    #Get total number of sequences in the current system
    n <- length(unlist(curr.otu_map));
    #Calculate shared mapping
    writeLines(paste0(date(), " => ", data[[sr]][[make.comparisons[i,1]]]$name, " vs ", data[[sr]][[make.comparisons[i,2]]]$name));
    writeLines(paste0(date(), " => Calculating current shared mapping..."));
    current.shared <- mclapply(curr.otu_map, function (map.2) {
      tmp <- tabulate(curr.seq_map[map.2]);
      out <- list();
      out$n.ij <- as.numeric(tmp[tmp != 0]);
      out$sizes.1 <- sizes.1[which(tmp != 0)];
      out
    }, mc.cores = PARAM$use.cores);
    #Compute Mutual Information & Entropies
    #=> equations (4) & (5) in Vinh et al., 2009
    writeLines(paste0(date(), " => Calculating Mutual Information and Entropies..."));
    mutual.information <- sum(unlist(mclapply(current.shared, function (shared) {sum((shared$n.ij / n) * log((shared$n.ij * n) / (sum(shared$n.ij) * shared$sizes.1)))}, mc.cores = PARAM$use.cores)));
    entropy.1 <- sum((sizes.1[sizes.1 != 0] / n) * log(sizes.1[sizes.1 != 0] / n));
    entropy.2 <- sum((sizes.2[sizes.2 != 0] / n) * log(sizes.2[sizes.2 != 0] / n));
    #Calculate Expected Mutual Information (EMI)
    writeLines(paste0(date(), " => Calculating Expected Mutual Information..."));
    tab.sizes.1 <- as.numeric(tabulate(sizes.1)); idx.sizes.1 <- which(tab.sizes.1 != 0);
    tab.sizes.2 <- as.numeric(tabulate(sizes.2)); idx.sizes.2 <- which(tab.sizes.2 != 0);
    expected.mutual.information <- sum(unlist(mclapply(
      idx.sizes.1,
      function (a.i) {
        sum(tab.sizes.1[a.i] * tab.sizes.2[idx.sizes.2] * unlist(lapply(
          idx.sizes.2,
          function (b.j) {
            n.ij <- as.numeric(seq(max(a.i + b.j - n, 1), min(a.i, b.j)));
            num.a.i <- as.numeric(a.i); num.b.j <- as.numeric(b.j);
            sum((n.ij / n) * log((n.ij * n) / (num.a.i * num.b.j)) * exp(lfactorial(num.a.i) + lfactorial(num.b.j) + lfactorial(n - num.a.i) + lfactorial(n - num.b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(num.a.i - n.ij) - lfactorial(num.b.j - n.ij) - lfactorial(n - num.a.i - num.b.j + n.ij)))
          }
        )))
      },
      mc.cores = PARAM$use.cores
    )));
    #Calcualte AMI
    res.ami.list[[sr]]$res[i] <- (mutual.information - expected.mutual.information) / (sqrt(entropy.1 * entropy.2) - expected.mutual.information);
    
    #Stop the timer
    writeLines(paste0(date(), "=> AMI ", res.ami.list[[sr]]$res[i]));
  }
}
#Export AMI results as tsv tables
write.table(res.ami.list[["v13"]]$res, file=paste0(PARAM$folder.results, "ami.v13.tsv"), quote=F, sep="\t", col.names=NA);
write.table(res.ami.list[["v35"]]$res, file=paste0(PARAM$folder.results, "ami.v35.tsv"), quote=F, sep="\t", col.names=NA);

################################################################################
################################################################################
q()












