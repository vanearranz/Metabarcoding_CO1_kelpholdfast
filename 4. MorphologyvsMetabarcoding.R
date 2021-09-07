########################################################
###### Script to compare Morphological taxa and Metabarcoding taxa 
########################################################

# This time with the new taxonomic assignments from a BLast setting without max_target_seqs 

require(VennDiagram)
require(dplyr)
require(gridExtra)
require(phyloseq)
require(grid)
require(eulerr)

### Metabarcoding taxa list ####

# To genererate OTU_nf_ass.out without controls and the bioinformatically combined synthetic samples, we have to run ASV_f1 with the function Edesign.merging (found in the Script 1.Edesign_functions_results.R)
OTU_nf_ass = subset_taxa(OTU_nofilter, Domain=="d__Eukaryota")
OTU_nf_ass.out  <- Edesign.merging(phyloseq.obj.raw = OTU_nf_ass, sample.data = sample_data_tomerge)

OTU_nf_ass.out

# Extracted only the 8 pooled samples from Matheson Bay (MTB05, MTB09, MTB11, MTB14...)

MTB_OTU_pool_samples <- subset_samples(OTU_nf_ass.out, sample_names(OTU_nf_ass.out)=="MTB05" | sample_names(OTU_nf_ass.out)=="MTB09" | sample_names(OTU_nf_ass.out)=="MTB11" |sample_names(OTU_nf_ass.out)=="MTB14"|sample_names(OTU_nf_ass.out)=="MTB20"|sample_names(OTU_nf_ass.out)=="MTB21"|sample_names(OTU_nf_ass.out)=="MTB27"|sample_names(OTU_nf_ass.out)=="MTB30")
sample_names(MTB_OTU_pool_samples)

# Prune taxa not present in the 8 samples 
MTB_OTUpool_samples_filtered = prune_taxa(!taxa_sums(MTB_OTU_pool_samples) < 1, MTB_OTU_pool_samples)
otu_table(MTB_OTUpool_samples_filtered)

write.csv(tax_table(MTB_OTUpool_samples_filtered), "MTBsamples_OTUassigned_taxtable.csv")

# Refine the metabarcoding assignmetns : 
## Remove Eukaryota only assigned to Eukaryota
## Remove the Domain column
## What should we do with the parasites? I removed them now : Paramoeba, Neoparamoeba, Cunia...
## Phaeophyceae assgned to Phylum Ochrophyta 
## Oomycetes is unaccepted Class from the Phylum Oomycota
## Dinophyceae is a class from the PHylum Myzozoa
## Pycnococcaceae is a Family from the Phylum Chlorphyta (reformatted the taxonomy according to WORMS) 
 
## AFTER MARINE CHECKER using - Worms_marine_checker.R
# Remove Phoridae as non-marine but leave Diptera because it can be marine (WORMS)

Metabarcoding_8MTB_newtaxonomy_clean_06082021 <- read.delim("Resources/Metabarcoding_8MTB_newtaxonomy_clean_06082021.txt")#, stringsAsFactors=TRUE)

## Unique rows and reformatted blanks
metab <- Metabarcoding_8MTB_newtaxonomy_clean_06082021
str(metab)
metab[metab==" "] <- NA
metab[metab==""] <- NA
str(metab)

# Metabarcoding unique taxa by phyllum results 

spp_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(distinct_species=n_distinct(Species, na.rm = TRUE))
genus_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(distinct_genus=n_distinct(Genus, na.rm = TRUE))
family_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(distinct_family=n_distinct(Family, na.rm = TRUE))
order_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(distinct_order=n_distinct(Order, na.rm = TRUE))
class_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(distinct_class=n_distinct(Class, na.rm = TRUE))
Phylum_distict_levels <- cbind(spp_per_phylum_metabarc, genus_per_phylum_metabarc$distinct_genus, family_per_phylum_metabarc$distinct_family, order_per_phylum_metabarc$distinct_order, class_per_phylum_metabarc$distinct_class)

write.csv(Phylum_distict_levels, "By_phyllum_METABARC.csv", row.names = FALSE)

# Metabarcoding ASSIGNED OTUs taxa by phyllum results

otu_spp_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(species=NROW(na.omit(Species)))
otu_genus_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(genus=NROW(na.omit(Genus)))
otu_family_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(family=NROW(na.omit(Family)))
otu_order_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(order=NROW(na.omit(Order)))
otu_class_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(class=NROW(na.omit(Class)))
otu_phylum_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(phylum=NROW(na.omit(Phylum)))

Phylum_distict_OTUlevels <- cbind(otu_spp_per_phylum_metabarc, otu_genus_per_phylum_metabarc$genus, otu_family_per_phylum_metabarc$family, otu_order_per_phylum_metabarc$order, otu_class_per_phylum_metabarc$class, otu_phylum_per_phylum_metabarc$phylum)

write.csv(Phylum_distict_OTUlevels, "OTU_By_phyllum_METABARC.csv", row.names = FALSE)

otu_spp_per_phylum_metabarc <- metab %>% group_by(Phylum) %>% summarize(species=NROW(Species))

### Morphological taxa list ####

# Merge Genus + species (when both are correct)

# Remove from Family column families with : Triticellidae?, indet. , ?
# Remove from Order column orders with : (subClass Thecatae),  (s-class)
# Crustacea is not a Class, so edit the proper Class according to the WORMS classification
# Amphipoda, Cumacea, Decapoda, Isopoda, Mysida (accepted name of Mysidacea) and Tanaidacea are Orders that belong to the Class Malacostraca. Crustacea is a sub-phyllum we don't include them here.
# Balanidae is a Family from the Order Sessilia and Class Hexanauplia (Change Crustacea as previously assigned to its Class)
# Copepoda is an infraclass. I set it up as Harpacticoida as it is the most common copedpoda found in metabarcoding and belong to the class Hexanauplia

Morphology_taxa <- read.csv("Resources/Morphology_9Leigh_clean.txt")

morph <- Morphology_taxa
morph[morph==" "] <- NA
morph[morph=="  "] <- NA
morph[morph==""] <- NA
str(morph)

# Morphology unique taxa by phyllum results

spp_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_species=n_distinct(Species.name, na.rm=TRUE))
genus_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_genus=n_distinct(Genus, na.rm = TRUE))
family_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_family=n_distinct(Family, na.rm = TRUE))
order_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_order=n_distinct(Order, na.rm = TRUE))
class_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_class=n_distinct(Class, na.rm = TRUE))
phylum_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(distinct_class=n_distinct(Phylum, na.rm = TRUE))

Phylum_distict_levels <- cbind(spp_per_phylum_morph, genus_per_phylum_morph$distinct_genus, family_per_phylum_morph$distinct_family, order_per_phylum_morph$distinct_order, class_per_phylum_morph$distinct_class)

write.csv(Phylum_distict_levels, "By_phyllum_morph.csv", row.names = FALSE)

# Morphology ASSIGNED OTUs taxa by phyllum results

otu_spp_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(species=NROW(na.omit(Species.name)))
otu_genus_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(genus=NROW(na.omit(Genus)))
otu_family_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(family=NROW(na.omit(Family)))
otu_order_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(order=NROW(na.omit(Order)))
otu_class_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(class=NROW(na.omit(Class)))
otu_phylum_per_phylum_morph <- morph %>% group_by(Phylum) %>% summarize(phylum=NROW(na.omit(Phylum)))

Phylum_distict_OTUlevels <- cbind(otu_spp_per_phylum_morph, otu_genus_per_phylum_morph$genus, otu_family_per_phylum_morph$family, otu_order_per_phylum_morph$order, otu_class_per_phylum_morph$class, otu_phylum_per_phylum_morph$phylum)

write.csv(Phylum_distict_OTUlevels, "OTU_By_phyllum_MORPH.csv", row.names = FALSE)

### Compare both db at SPECIES level #### 

# Remove from column with : indet. , sp. sp1.. 
morph_clean <- Morphology_taxa
morph_clean[morph_clean==" "] <- NA
morph_clean[morph_clean=="  "] <- NA
morph_clean[morph_clean==""] <- NA
str(morph)

metab_sp <- unique(na.omit(metab$Species))
morph_sp <- unique(na.omit(morph_clean$Species.name))

spp_both <- metab_sp[metab_sp %in% morph_sp]
spp_only_metab <- metab_sp[!metab_sp %in% morph_sp]
spp_only_morph <- morph_sp[!morph_sp %in% metab_sp]
length(spp_both)
length(spp_only_metab)
length(spp_only_morph)

## Compare GENUS level 
metab_genus <- unique(na.omit(metab$Genus))
morph_genus <- unique(na.omit(morph_clean$Genus))

genus_both <- metab_genus[metab_genus %in% morph_genus]
genus_only_metab <- metab_genus[!metab_genus %in% morph_genus]
genus_only_morph <- morph_genus[!morph_genus %in% metab_genus]
length(genus_both)
length(genus_only_metab)
length(genus_only_morph)

## Compare FAMILY level 
metab_family <- unique(na.omit(metab$Family))
morph_family <- unique(na.omit(morph_clean$Family))

family_both <- metab_family[metab_family %in% morph_family]
family_only_metab <- metab_family[!metab_family %in% morph_family]
family_only_morph <- morph_family[!morph_family %in% metab_family]
length(family_both)
length(family_only_metab)
length(family_only_morph)

## Compare ORDER level 
metab_order <- unique(na.omit(metab$Order))
morph_order <- unique(na.omit(morph_clean$Order))

order_both <- metab_order[metab_order %in% morph_order]
order_only_metab <- metab_order[!metab_order %in% morph_order]
order_only_morph <- morph_order[!morph_order %in% metab_order]
length(order_both)
length(order_only_metab)
length(order_only_morph)

## Compare CLASS level 
metab_class <- unique(na.omit(metab$Class))
morph_class <- unique(na.omit(morph_clean$Class))

class_both <- metab_class[metab_class %in% morph_class]
class_only_metab <- metab_class[!metab_class %in% morph_class]
class_only_morph <- morph_class[!morph_class %in% metab_class]
length(class_both)
length(class_only_metab)
length(class_only_morph)

## Compare PHYLLUM level 
metab_phylum <- unique(na.omit(metab$Phylum))
morph_phylum <- unique(na.omit(morph_clean$Phylum))

phylum_both <- metab_phylum[metab_phylum %in% morph_phylum]
phylum_only_metab <- metab_phylum[!metab_phylum %in% morph_phylum]
phylum_only_morph <- morph_phylum[!morph_phylum %in% metab_phylum]
length(phylum_both)
length(phylum_only_metab)
length(phylum_only_morph)

## Few Venn diagrams

Species_venn <- plot(euler(c("A"=63,"B"=107,"A&B"=2)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))
Genus_venn <- plot(euler(c("A"=89,"B"=115,"A&B"=6)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))
Family_venn <- plot(euler(c("A"=72,"B"=105,"A&B"=19)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))
Order_venn  <- plot(euler(c("A"=53,"B"=31,"A&B"=15)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))
Class_venn <- plot(euler(c("A"=22,"B"=7,"A&B"=13)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))
Phylum_venn <-plot(euler(c("A"=7,"B"=3,"A&B"=11)),quantities = TRUE,fills=c("#6B8AB6","#EEEEB1","#B6D6AF"))

gridExtra::grid.arrange(Species_venn, Genus_venn,Family_venn, Order_venn, Class_venn, Phylum_venn, ncol=3, nrow=2)


##### EXTRA  ###########################################################################
################################################################################
################################################################################
### Reviewer 1 mentioned that we can add some numbers for comparing the visual surveys species and what sequences are available in the reference database for those species 

file<-"Resources/seqnames_MARES_June2021_nobarcode.txt" ## File containing sequence names
mares_names<-readLines(file)#,sep=" ",stringsAsFactors=FALSE,head=F)
mares_names<-gsub("^ *|(?<= ) | *$", "", mares_names, perl = TRUE)
mares_names<-strsplit(mares_names," ")
mares_names<-lapply(mares_names,head,3)
mares_names<-as.data.frame(do.call("rbind",mares_names))

library(stringr)
convert2name<-function(seqnames){
  seqnames<-paste(seqnames$V2,seqnames$V3)
  seqnames<-trimws(gsub("\\w*[0-9]+\\w*\\s*", "", seqnames))
  seqnames<-trimws(gsub("\\w*\\.+\\w*\\s*", "", seqnames))
  seqnames
}
mares_seq_names <-convert2name(mares_names)
mares_uniseq_names <-unique(mares_seq_names)

# Compare both MARES and Morphology 

# At species level 
intersect(mares_uniseq_names, morph_sp)

### 17 SPECIES IN COMMON
#[1] "Evechinus chloroticus"    "Lysidice ninetta"         "Onithochiton neglectus"   "Caprella equilibra"      
#[5] "Ophiothrix aristulata"    "Patiriella regularis"     "Hiatella arctica"         "Plumularia setacea"      
#[9] "Ophiopteris antipodum"    "Beania plurispinosa"      "Arachnopusia unicornis"   "Beania magellanica"      
#[13] "Calloporina angustipora"  "Crella incrustans"        "Cookia sulcata"           "Crepidacantha crinispina"
#[17] "Galeopsis porcellanicus

#Compare at genus level 
genus_foundinmares <- intersect(mares_uniseq_names, morph_genus)
length(intersect(mares_uniseq_names, morph_genus_names))
#40
 
#[1] "Hiatella"     "Platynereis"  "Caprella"     "Callyspongia" "Lamellaria"   "Pagurus"      "Nereis"       "Ophiactis"   
#[9] "Alpheus"      "Didemnum"     "Barbatia"     "Pilumnus"     "Plumularia"   "Cnemidocarpa" "Balanus"      "Terebellides"
#[17] "Flabelligera" "Photis"       "Eatonina"     "Ophiothrix"   "Dysidea"      "Iphimedia"    "Branchiomma"  "Eunice"      
#[25] "Lysidice"     "Orthopyxis"   "Leucothoe"    "Orchomene"    "Trochus"      "Emarginula"   "Scissurella"  "Eatoniella"  
#[33] "Epitonium"    "Polycheria"   "Melita"       "Crella"       "Microporella" "Fenestrulina" "Lepidonotus" 


## Compare MARES with all species found in the holdfast with METABARCODING 

metab_species_names <-"./metab_kelp_holdfast_allspecies.txt"
metab_species_names <- readLines(metab_species_names)
metab_species_names <- unique(metab_species_names)

metabspecies_foundinmorph <- intersect(morph_names, metab_species_names)

# 5 species 
#[1] "Beania magellanica"      "Beania plurispinosa"     "Calloporina angustipora" "Evechinus chloroticus"   "Plumularia setacea"  