# Metabarcoding_CO1_kelpholdfast
This workflow is specific to analysing eukaryote diversity from cDNA metabarcoding of kelp holdfasts using CO1 region from raw reads to biostatiscal analysis included in the manuscript "Metabarcoding hyperdiverse marine communities in temperate kelp forests: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre. 

CO1 region was amplified using mlCOIintF-XT (refs) and jgHCO2198 (refs) primers. The product was aproximately 313bp obtained by Ilumina paired-end sequencing. 

Sample preparation, DNA extraction and PCR amplification related to this work can be found here (add link).

The scripts are designed to be run using a Linux OS, and were developed on Ubuntu 16.04. 

## Requirements : 

- The contents of this repo
- Raw sequence data : CHECK WHICH ONE IS THE EASIEST WAY TO DOWNLOAD IT READY TO USE FOR THIS PIPELINE. Now I have them in a local folder name processed_trimmed. 
All FASTQ sequence files are available from the National Center for Biotechnology Information short-read archive database (Bioproject: PRJNA638997, Biosamples: SAMN15220525-SAMN15220620).
- QIIME2 version 2019.4 https://docs.qiime2.org/2020.8/install/ 
- Biom http://biom-format.org/
- MARES reference sequences database : https://osf.io/4f8mk/ 
- BLASTn  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- MEGAN6 -  Metagenome Analyzer https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/megan6/
- VSEARCH https://github.com/torognes/vsearch
- R https://www.r-project.org/


## Getiing started

### Activate QIIME2 

**Citation:** Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). https://doi.org/10.1038/s41587-019-0209-9

```
conda activate qiime2-2019.4
``` 
### Import raw sequence data
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path processed_trimmed/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza


# paired-end-demux.qza is an artifact file '.qza' that has the information stored for your raw fastq sequences. We can summarize and visualize the output by transforming to '.qzv'. Check the output : https://view.qiime2.org/

qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv

qiime tools view demux-paired-end.qzv
```
## Remove primers : cutadapt

We use [cutadapt](https://github.com/qiime2/q2-cutadapt) to remove the primers

**Citation:** Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1):pp–10, 2011. https://doi:10.14806/ej.17.1.200.

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza \
  --p-front-f GGWACWRGWTGRACWNTNTAYCCYCC \
  --p-front-r TANACYTCNGGRTGNCCRAARAAYCA \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqs.qza \
  --verbose
  
qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv
  
qiime tools view trimmed-seqs.qzv
```

## Denoise, chimera removal and clustering in ASVs: DADA2

[DADA2](https://github.com/qiime2/q2-dada2) : Pair-end joining, dereplication, chimera filtering and clustering in ASVs 

**Citation:** Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869

We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).

We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.

ASVs are inferred by a *de novo* process in which biological sequences are discriminated from errors on the basis of, in part, the expectation that biological sequences are more likely to be repeatedly observed than are error-containing sequences. As a result, ASV inference cannot be performed independently on each read—the smallest unit of data from which ASVs can be inferred is a sample. However, unlike *de novo* OTUs, ASVs are consistent labels because ASVs represent a biological reality that exists outside of the data being analyzed: the DNA sequence of the assayed organism. Thus, ASVs inferred independently from different studies or different samples can be validly compared.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --p-n-threads 4 \
  –-p-chimera-method consensus \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza
  
qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv 
qiime tools view table-dada2.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv \
qiime tools view rep-seqs-dada2.qza

# From this visualization we can download the FASTA FILE that we will use for taxonomic assignment (see Taxonomic assignment section)

qiime metadata tabulate \
--m-input-file denoising-stats-dada2.qza \
--o-visualization denoising-stats-dada2.qzv
```

### Export ASV table and representative sequences 

```
# Example ASV table without filtering 
qiime tools export
-- input-path table-dada2.qza \
--output-path ASVtable/

#Export representative sequences
qiime tools export  \
--input-path rep-seqs-dada2.qza \
--output-path rep-seq-ASV.fasta
```

Convert the ASV table to tsv in Biom 

```
biom convert
-i ASVtable/feature-table.biom
-o ASVtable/ASV-frequency-table.tsv  --to-tsv
```

## Taxonomic assignment 

ASVs passing the quality control and filtering thresholds (rep-seq-ASV.fasta) were taxonomically assigned using the MARES reference sequence database. 

[MARES](https://www.nature.com/articles/s41597-020-0549-9) is the most comprehensive CO1 reference database for marine eukaryotes available, and provides users the ability to retain taxa that cannot be assigned at the species level, but can be assigned at higher taxonomic levels. 

**Citation:** Arranz, V., Pearman, W.S., Aguirre, J.D. et al. MARES, a replicable pipeline and curated reference database for marine eukaryote metabarcoding. Sci Data 7, 209 (2020). https://doi.org/10.1038/s41597-020-0549-9


**To use MARES reference database:**

Download MARES_NOBAR from https://osf.io/4f8mk/
*Note : You can also used you MARES_BAR.*

Place MARES_NOBAR_BOLD_NCBI_sl_reformatted.fasta in the main folder or include the path after -db 

### Build the database 

We first create a blast database from MARES reference sequence database
```
makeblastdb -in MARES_NOBAR_BOLD_NCBI_sl_reformatted.fasta -dbtype nucl -parse_seqids
```

### Blast each sequence against MARES reference sequences database : Blastn 

Then, we performed a BLASTn against MARES reference database with an e-value of 1-60 for high-quality matches and max_target_seqs equal to 10. 

```
blastn -db MARES_NOBAR_BOLD_NCBI_sl_reformatted.fasta -query rep-seq-ASV.fasta -evalue 1e-60 -max_target_seqs 10 -outfmt 5 -out MARES_MEGAN.txt -num_threads 12

```

### Use LCA algorithm for taxonomic assignment : MEGAN6

We used MEGAN 6.18.9 for taxonomic assignment within the NCBI taxonomy framework using the default Lowest Common Ancestor (LCA) algorithm parameters. 

**Citation** : Huson DH, Beier S, Flade I, Górska A, El-Hadidi M, et al. (2016) MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLOS Computational Biology 12(6): e1004957. https://doi.org/10.1371/journal.pcbi.1004957

Launch MEGAN 6.18.9

Import the blast output and the fasta file used as the blast query into MEGAN (File → Import From BLAST) : MARES_MEGAN.txt and rep-seq-ASV.fasta 

Apply the following LCA settings:

```
min score 98 
max expected 0.00000001 
min % ID 75 
top % 10 
min support % 0 (off) 
min support 1 
```

Select level of taxonomy to view, possibly use multiple different levels, e.g. species, genus, family
File -> Export -> Text (csv) Format
Choose: readName_to_taxonPathKPCOFGS

-	It has a percent value at the end of the line. This refers to the percentage of high scoring alignments for the given read that map to the last taxon on the path. It has nothing to do with the percentage used in the weighted LCA.
-	It only reports taxa in the path that have an official KPCOFGS rank. Intermediate nodes that have no taxonomic rank, or one that does not belong to KPCOFGS, are suppressed
-	Each node is prefixed by letter__ to indicate the rank, e.g. g__ for genus, s__ for species : we can edit this later to have only the name

Save into the main folder as assigned_seqs-MARES-ex.txt 


**To create the OTU_taxonomy file to use in downstream analysis (see Statistical analysis):**
- Export it as ReadName_to_taxonPathKPCOFGS / ASSIGNED
- Import it into excel and remove the columns with the percentage *No need in the last version of MEGAN 
- The unknowns were giving problems. Edited the excel to get rid of them and replace them for the taxon known 
- The ASVs not assigned have to be set as 'd_unassigned' because is not working otherwise
- Generate the names with a semicolon ; separation
- Create OTU_ID in the first column and the TaxonPAth in the other one
- Save as .csv or .txt -> **taxonomy_edited_8ranks.txt**

## Create ASV/OTU Phyloseq objects 

We import the ASV table in R as [phyloseq](https://joey711.github.io/phyloseq/) object with the taxonomy associated, sample metadata and reference sequences for posterior statistical analysis. 

```
library("phyloseq")

# ASV table without filtering 

ASV_table <- read_delim("ASVtable.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

## Add the taxonomy 

tax_table_new_edited_8ranks <- read.delim("/Volumes/NO NAME/eDNA revisions/tax_table_new_edited_8ranks.txt", row.names=1)
str(tax_table_new_edited_8ranks)

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}
newtaxonomy <-matrix.please(tax_table_newtaxonomy_8ranks)

tax_table_newtaxonomy_8ranks <- tax_table(newtaxonomy)


# Add the reference sequences
reference_seqs0 <- readDNAStringSet(file = "rep-sequences-nofilt.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
REFS00 = refseq(reference_seqs0)

# Add the sample data file 
sample_data_batches <- read.csv("~/Desktop/Metabarcoding_CO1_kelpholdfast/sample_data_11varX96samples.csv")
sampledata00 = sample_data(data.frame(sample_data_batches, row.names = sample_names(ps_asv00)))
str(sampledata00)

# Create the phyloseq object 
ASV <- phyloseq(otu_table(ASV_nofilter), sample_data(sampledata00), refseq(reference_seqs0), tax_table(tax_table_newtaxonomy_8ranks))

```

## Extract possible contaminants and tag switching normalisation 

We used decontam package in R to identify and extract blank contaminants using the negative controls. 

**Citation:**  Davis NM, Proctor D, Holmes SP, Relman DA, Callahan BJ (2017). “Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data.” bioRxiv, 221499. https://doi.org/10.1101/221499 


```
library("phyloseq")
library("decontam")

#Remove Seawater control 
ASV_nf.noSC <-  subset_samples(ASV, Sample != "SC")

# Identify and extract blank contaminants using the negative controls. The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants ###

sample_data(ASV_nf.noSC)$is.neg <- sample_data(ASV_nf.noSC)$Sample_or_Control == "Control"
contam.comb.all <- isContaminant(ASV_nf.noSC, method="combined", neg="is.neg", threshold=0.5, conc="Concentration")

# Extract contaminants
ASV_nf.noSC.nocon <- prune_taxa(!contam.comb.all$contaminant, ASV_nf.noSC)

# Extract negative control samples 
ASV_nf.noSC.nocon.nc <- subset_samples(ASV_nf.noSC.nocon, Sample_or_Control=="True Sample")

write.csv(otu_table(ASV_nf.noSC.nocon.nc),file = "ASVtable_nf_nosc_noc_nocon.csv")

```

Tag switching correction. We used the R Script included in Resources/owi_renormalize.R ADD THE SCRIPT THERE!!!! (https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/owi_renormalise.R) It sort the samples by abundance of each ASV and eliminate the reads of the samples corresponding to a cumulative frequency of less than 3% for each particular ASV . 

**Citation:**  Wangensteen OS, Turon X (2016) Metabarcoding techniques for assessing biodiversity of marine animal forests. Marine Animal Forests. The Ecology of Benthic Biodiversity Hotspots, eds Rossi S, Bramanti L, Gori A, Orejas C (Springer International Publishing).

```
RScript owi_renormalize.R -i ASVtable_nf_nosc_noc_nocon.tsv -o ASVtable_tsc.tsv -c 0.97 -s 2 -e 88
```

Rename samples names from . to - in the output table ASVtable_tsc.tsv. In R : 

```
#After tag switching normalization
ASVtable_tsc <- read.delim("~/ASVtable_tsc.tsv",sep = "\t", header=TRUE,as.is=TRUE, row.names = 1, check.names = FALSE)
# delete the last 3 columns of total counts 
ASVtable_tsc <- ASVtable_2lulu2[-c(88:90)]
# Remove rows that sum columns are 0
ASVtable_tsc <- ASVtable_tsc[as.logical(rowSums(ASVtable_tsc != 0)), ]


ASV_nofilter <- otu_table(ASVtable_tsc_all, taxa_are_rows = TRUE)

# Create Phyloseq object with ASV table after tag switching normalisation 
ASV_nofilter <- phyloseq(otu_table(ASV_nofilter), sample_data(sampledata00), refseq(reference_seqs0), tax_table(tax_table_newtaxonomy_8ranks))

```


## Clustering ASVs into OTUs : VSEARCH 

We used Vsearch to cluster the ASVs into Operational Taxonomic Units (OTUs). We chose 97% of similarity for CO1 mitochondrial region.

**Citation:**  Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584

```
vsearch --cluster_size rep-sequences-ASV.fasta  --id 0.97 --uc clustering-results.uc -msaout sequences-outs
```

Save the clustering results.uc as .csv

I edited manually deleting the consensus rows C (because they are the same as the S of the clusters) -> clustering-results.csv

Change the header names as ASV_ID to merge in R with the previously taxonomy assigned to ASVs.

By running the following script in R, we created the OTU table by combining the previous ASV table and Vsearch clustering results by ASV_ID column and then.  
 
``` {r}
clustering.results <- read.csv("~/Desktop/Edesign manuscript/eDNA revisions/clustering-results.csv")
ASV_OTU_table <- merge(ASVtable_tsc_all_clean, clustering.results, by.x="ASV_ID", by.y="ASV_ID")
write.csv(ASV_OTU_table, file = "ASV_OTU_table.csv")
```

Manually edited again the ASV_OTU_table.csv to create the ASV_OTU_tabletocollapse.csv
- Sort by H and  S : copy all the ASV column of the S into the OTU column (they are the consensus of the clusters)
- Sort by Cluster number
- Delete all the columns from vsearch : taxonomy and ASV_ID

```
OTU_nofilter_tocollapse <-  read.csv("~/Desktop/Edesign manuscript/eDNA revisions/ASV_OTU_tabletocollapse.csv", check.names = FALSE)
OTU_nofilter_tocollapse$OUT97_ID <- as.factor(OTU_nofilter_tocollapse$OUT97_ID)
#Collapse the OTUs with the same name (all ASVs that were collapse into OTUs)
OTU_nofilter_collapsed <- OTU_nofilter_tocollapse %>% group_by(OUT97_ID) %>% summarise_all(funs(sum))

OTU_nofilter_collapsed$OUT97_ID <- as.character(OTU_nofilter_collapsed$OUT97_ID)
OTU_nofilter_collapsed <- as.data.frame(OTU_nofilter_collapsed)
rownames(OTU_nofilter_collapsed) <- OTU_nofilter_collapsed[,1]
OTUtable_tsc <- OTU_nofilter_collapsed[,-1]


# Create Phyloseq object with OTU table after tag switching normalisation 
OTU_nofilter <- otu_table(OTUtable_tsc, taxa_are_rows = TRUE)

OTU_nofilter <- phyloseq(otu_table(OTU_nofilter), sample_data(sampledata00), refseq(reference_seqs0), tax_table(tax_table_newtaxonomy_8ranks))

```

## Refining the datasets for downstream analysis 

At this stage we perform different strategies to remove sequencing errors and artifacts to explore the effects of a number of different filtering thresholds (site-occupancy vs. percentahe cut-off) in the biodiversity estimates. We also included a non-filtered table for our downstream analysis.


### LULU 

ASV/OTU table curation combining similarity and co-occurence patterns. 

**Citation:** Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature communications, 8(1), 1-11.

```
###  LULU : ASV/OTU table curation based on similarity and co-occurence rates (Froslev et al. 2017)
library("magrittr")
library("lulu")

# I need : 
######## a. ASV/OTU table 

# ASV table
ASVtable_tsc 
#OTU table 
OTUtable_tsc
####### b. Produce a match list from the fasta file with the sequences

# Extract the ref_seqs from the phyloseq object
write.csv(refseq(ASV_nf.noSC.nocon.nc), ref_seqs.csv)
#I convert csv to fasta in a website - rep_sequences-ASV.fasta
```

In the TERMINAL with BLASTN. First produce a blastdatabase with the ASV/OTUs reference sequences 

```
makeblastdb -in rep-sequences-ASV.fasta -parse_seqids -dbtype nucl
```
Then blast the OTUs against the database
```
blastn -db rep-sequences-ASV.fasta -outfmt '6 qseqid sseqid pident' -out match_list2.txt -qcov_hsp_perc 80 -perc_identity 84 -query rep-sequences-ASV.fasta
```

Back in R 

```
matchlist2 <- read.table("match_list2.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

### Run LULU to obtained the new OTU table
##### WITH tag switching normalization 
lulu_curated_result_ASV_tsc <-lulu(ASVtable_tsc, matchlist2)
#Total of 7019 (excluding rows that sum columns are 0)
lulu_curated_result_ASV_tsc$curated_count
#4577
lulu_curated_result_ASV_tsc$discarded_count
#2442

# Create a Phyloseq object with the new OTU table with LULU 
ASVtable_lulu <- as.matrix(lulu_curated_result_ASV_tsc$curated_table)
str(ASVtable_lulu)

ASV_lulu <- otu_table(ASVtable_lulu, taxa_are_rows = TRUE)
ASV_lulu <- phyloseq(otu_table(ASV_lulu), sample_data(ASV_nf), refseq(ASV_nf), tax_table(ASV_nf))
```

### Minimum read abundance filtering 

One of the most broadly employed strategies to remove artefactual sequences is to discard sequences with copy numbers under a certain threshold.  We can do many diferent filters at this stage, we used -p-min-frequency to remove low abundance features of less than 0.003%, 0.01% and 0.05% across all samples.  

```
ADD THE SCRIPT FOR THIS 
```


## Statistical analysis 

Merge samples and perform linear models 

The remaining Phyloseq objects are included in Resources/[Edesign.Rdata](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/Edesign.Rdata) file ready to perform next step.  

- Edesign.cleaning : to merge samples to generate the bioinformatically merge samples. 

We used linear models implemented in R v3.6.1 (R Core Team, 2019) to assess how our laboratory and bioinformatic decisions influenced biodiversity estimates for the holdfast communities. 

We created 4 different ASV tables by applying 3 abundance filtering thresholds and without filtering : ASVnf, ASV39, ASV130 and ASV652.
Each of the ASV tables was converted into OTU by clustering at 97% similarity the ASVs : OTUnf, OTU39, OTU130 and OTU652. 

- Edesign.LM.pipeline : uses the contrast matrix we used for our tests and run the linear models with the rarefied presence/absence matrix to assess how our laboratory and bioinformatic decisions influenced biodiversity estimates for the holdfast communities. 

We run the custom R script with the Phyloseq objects created in the previous steps and contained in [Edesign.Rdata](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/Edesign.Rdata) file. 

The results of the linear models were included in Figure 2 and Supporting information Table S2 of the manuscript "Metabarcoding hyperdiverse marine communities in temperate kelp forests: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre. 





