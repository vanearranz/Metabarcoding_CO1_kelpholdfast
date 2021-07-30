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
- 

## Create ASV/OTU Phyloseq objects 

We used a combination of Biom (in the terminal) and R Scripts, to import each ASV/OTU table in R as [phyloseq](https://joey711.github.io/phyloseq/) objects with the taxonomy associated, sample metadata and reference sequences for posterior statistical analysis. 

In R, combine ASV_table no filtered and taxonomy associated:
```
# Example with only ASV table without filtering 

asv_table_nf <- read_delim("ASVtable.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
taxonomy <- read_delim("taxonomy_edited_8ranks.txt",  "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
asvNF_taxonomy <- merge(asv_table_nf, taxonomy, by.x ="OTU_ID", by.y ="OTU_ID")
write.table(asvNF_taxonomy, "ASVnf_taxonomy.txt", sep="\t", row.names=FALSE) 
```

In the terminal, convert  the .txt into [Biom format](https://biom-format.org/) with the taxonomy as metadata 
```
biom convert -i ASVnf_taxonomy.txt -o asv.table.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy
```

Back in R, import the Biom table with taxonomy attached and create the phyloseq object
```
library("biomformat")
library("Biostrings")

x0 =read_biom("asv.table.biom")
#Add ASV table
asvmat00 = as(biom_data(x0), "matrix")
ASV00 = otu_table(asvmat00, taxa_are_rows=TRUE)
#Add taxonomy 
taxmat00 = as.matrix(observation_metadata(x0), rownames.force=TRUE)
TAX00 = tax_table(taxmat00)
# Add the reference sequences
reference_seqs0 <- readDNAStringSet(file = "rep-sequences-nofilt.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
REFS00 = refseq(reference_seqs0)

# Create the phyloseq object
ps_asv00 = phyloseq(OTU00, TAX00, REFS00)

# Add the sample data file 
sample_data_batches <- read.csv("~/Desktop/Metabarcoding_CO1_kelpholdfast/sample_data_11varX96samples.csv")
sampledata00 = sample_data(data.frame(sample_data_batches, row.names = sample_names(ps_asv00)))
str(sampledata00)

# Merge with phyloseq object 
ASV_nf <- merge_phyloseq(ps_asv00, sampledata00)
```

### Extract possible contaminants and tag switching normalisation 

- ADD R Script for this part!!! 

## Clustering ASVs into OTUs : VSEARCH 

We used Vsearch to cluster the ASVs into Operational Taxonomic Units (OTUs). We chose 97% of similarity for CO1 mitochondrial region.

**Citation:**  Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584

```
vsearch --cluster_size rep-seq-ASV.fasta  --id 0.97 --uc clustering-results.uc -msaout sequences-outs
```

Save the clustering results.uc as .csv

I edited manually deleting the consensus rows C (because they are the same as the S of the clusters) -> 97clusters-ASVintoOTUS.csv

Change the header names as ASV_ID to merge in R with the previously taxonomy assigned to ASVs.

By running the following script in R, we created the OTU table by combining the ASV table and Vsearch clustering results (97clusters-ASVintoOTUS.csv) by ASV_ID column and then, we attached the taxonomy.  
 
``` {r}
ASV_table <- read.delim("~/Desktop/Metabarcoding_CO1_kelpholdfast/ASV_table.csv")
97clusters.ASVintoOTUS <- read.delim("~/Desktop/Metabarcoding_CO1_kelpholdfast/97clusters-ASVintoOTUS.csv")
ASV_OTU_table <- merge(ASV_table, 97clusters.ASVintoOTUS, by.x="ASV_ID", by.y="ASV_ID")

taxonomy_edited_8ranks <- read.delim("~/Desktop/Metabarcoding_CO1_kelpholdfast/taxonomy_edited_8ranks.txt")
ASV_OTU_taxonomy_table <- merge(ASV_OTU_table, taxonomy_edited_8ranks, by.x="ASV_ID", by.y="OTU_ID")
write.csv(ASV_OTU_taxonomy_table, file = "ASV_OTU_taxonomy_table.csv")
```

Manually edited again the ASV_OTU_taxonomy_table.csv to create the OTU97_table.csv
- Sort by H and  S : copy all the ASV column of the S into the OTU column (they are the consensus of the clusters)
- Sort by Cluster number
- Delete all the columns from vsearch : taxonomy and ASV_ID

## Statistical analysis 

The remaining Phyloseq objects are included in Resources/[Edesign.Rdata](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/Edesign.Rdata) file ready to perform next step.  

- Edesign.cleaning : to merge samples to generate the bioinformatically merge samples. 

We used linear models implemented in R v3.6.1 (R Core Team, 2019) to assess how our laboratory and bioinformatic decisions influenced biodiversity estimates for the holdfast communities. 

We created 4 different ASV tables by applying 3 abundance filtering thresholds and without filtering : ASVnf, ASV39, ASV130 and ASV652.
Each of the ASV tables was converted into OTU by clustering at 97% similarity the ASVs : OTUnf, OTU39, OTU130 and OTU652. 

- Edesign.LM.pipeline : uses the contrast matrix we used for our tests and run the linear models with the rarefied presence/absence matrix to assess how our laboratory and bioinformatic decisions influenced biodiversity estimates for the holdfast communities. 

We run the custom R script with the Phyloseq objects created in the previous steps and contained in [Edesign.Rdata](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/Edesign.Rdata) file. 

The results of the linear models were included in Figure 2 and Supporting information Table S2 of the manuscript "Metabarcoding hyperdiverse marine communities in temperate kelp forests: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre. 

#merge samples and perform linear models 

## Abundance filtering 

We can do many diferent filters at this stage, we used -p-min-frequency to remove low abundance features of less than 0.003%, 0.01% and 0.05% across all samples, corresponding to 39, 130 and 651 sequence reads in our case. We also included a non-filtered table for our posterior analysis. 

```
# Example with 0.003% minimum number of sequences of total abundance (= 39 sequence reads)
qiime feature-table filter-features \
  --i-table table-dada2.qza  \
  --p-min-frequency 39 \
  --o-filtered-table filtered39-table-dada2.qza

qiime feature-table filter-seqs
--i-data rep-seqs-dada2.qza
--i-table filtered39-table-dada2.qza
--o-filtered-data filtered39-seqs-dada2.qza
```

