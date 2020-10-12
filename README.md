# Metabarcoding_CO1_kelpholdfast
This workflow is specific to analysing eukaryote diversity from cDNA metabarcoding of kelp holdfasts using CO1 region from raw reads to biostatiscal analysis included in the manuscript "Metabarcoding hyperdiverse marine communities in temperate kelp forests: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre. 

CO1 region was amplified using mlCOIintF-XT (refs) and jgHCO2198 (refs) primers. The product was aproximately 313bp obtained by Ilumina paired-end sequencing. 

Sample preparation, DNA extraction and PCR amplification related to this work can be found here (add link).

The scripts are designed to be run using a Linux OS, and were developed on Ubuntu 16.04. 

## Requirements : 

- This repo contents
- Raw sequence data : CHECK WHICH ONE IS THE EASIEST WAY TO DOWNLOAD IT READY TO USE FOR THIS PIPELINE. Now I have them in a local folder name processed_trimmed. 
All FASTQ sequence files are available from the National Center for Biotechnology Information short-read archive database (Bioproject: PRJNA638997, Biosamples: SAMN15220525-SAMN15220620).
- QIIME2 version 2019.4 https://docs.qiime2.org/2020.8/install/ 
- Biom (install instructions)
- MARES reference sequences database (intructions)
- BLASTn  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- MEGAN 
- VSEARCH (install instructions)
- R https://www.r-project.org/


## Getiing started

# Activate QIIME2 
```
conda activate qiime2-2019.4
``` 
# Import raw sequence data
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

DADA2 : Pair-end joining, dereplication, chimera filtering and clustering in ASVs 

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

### Export ASV table 

```
qiime tools export
-- input-path filtered25-table-dada2.qza
--output-path ASV25table/

qiime tools export  \
--input-path filtered25-seqs-dada2.qza \
--output-path Dada2-filtered25-repseqs
```

Convert the ASV table in Biom 

```
biom convert
-i ASV25table/feature-table.biom
-o ASV25table/ASV25-frequency-table.tsv  --to-tsv
```
## Taxonomic assignment 

ASVs passing the quality control and filtering thresholds were taxonomically assigned using the MARES_COI_NOBAR reference sequence database (Arranz, Pearman, Aguirre, & Liggins, 2019). 

MARES is the most comprehensive COI reference database for marine eukaryotes available, and provides users the ability to retain taxa that cannot be assigned at the species level, but can be assigned at higher taxonomic levels (Arranz, Pearman, Aguirre, & Liggins, 2020). 

We first performed a BLASTn (Altschul, Gish, Miller, Myers, & Lipman, 1990) with an e-value of 1-60 for high-quality matches and max_target_seqs equal to 1. Then, we used MEGAN 6.18.3 (Huson et al., 2016), for taxonomic assignment within the NCBI taxonomy framework using the default Lowest Common Ancestor algorithm parameters

# Assign each sequence to a taxon

BLASTN Prepare to mEGAN.SH 

# Annotate sequences with species

Launch MEGAN, import the blast output and the fasta file used as the blast query
Apply the following LCA settings:
min score 100
max expected 0.00000001
min % ID 97
top % 10
min support % 0 (off)
min support 1

Select level of taxonomy to view, possibly use multiple different levels, e.g. species, genus, family
File -> Export csv
Choose: readName_to_taxonName
Save into the folder 05_annotated
