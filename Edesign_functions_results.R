###### NEW RESULTS 2/08/21 ######
library(phyloseq)
library(decontam)
library(lme4)
library(lmerTest)
library(vegan)

############# FUNCTION TO MERGE SAMPLES AND CREATE THE BIOINFORMATICALLY COMBINED SAMPLES ######### 

Edesign.merging <- function(phyloseq.obj.raw, sample.data){
  
  ########## #### Extract samples for the Experimental design  
  phyloseq.obj.Edesign <- subset_samples(phyloseq.obj.raw, Edesign=="yes")
  
  # Create a phyloseq object from the otu_table, sample data that allow us to merge by factor (sample_data_tomerge)
  phyloseq.obj.norare <- phyloseq(otu_table(phyloseq.obj.Edesign), sample_data(data.frame(sample.data, row.names = sample_names(phyloseq.obj.Edesign))), tax_table(phyloseq.obj.Edesign), refseq(phyloseq.obj.Edesign)) 

  ########## Merging samples that will be compared in our experimental design
  
  ## Large + Small fractions (subset BCF, merge samples by Sample and merge with the previous phyloseq) 
  Lrg.Sml <- subset_samples(phyloseq.obj.norare, BCF=="BCF")
  BCF <- merge_samples(Lrg.Sml,group = "Smp", fun = sum)
  sample_names(BCF) <- c("MTB09_BCF", "MTB11_BCF", "MTB21_BCF")
  
  phyloseq.obj_norare.BCF <- merge_phyloseq(phyloseq.obj.norare, BCF)
  
  ## Ea + Eb (subset BCE, merge samples by Sample and merge with the previous phyloseq) 
  Ea.Eb <- subset_samples(phyloseq.obj.norare, BCE=="BCE")
  BCE <- merge_samples(Ea.Eb,group = "Smp", fun = sum)
  sample_names(BCE) <- c("MTB05_BCE", "MTB09_BCE", "MTB11_BCE", "MTB14_BCE", "MTB20_BCE", "MTB21_BCE", "MTB27_BCE", "MTB30_BCE")
  
  phyloseq.obj_norare.BCF.BCE <-merge_phyloseq(phyloseq.obj_norare.BCF, BCE)
  
  ## Px+Py (subset BCP, merge samples by Sample and merge with the previous phyloseq)
  Px.Py <- subset_samples(phyloseq.obj.norare, BCP=="BCP")
  BCP <- merge_samples(Px.Py,group = "SF", fun = sum)
  sample_names(BCP) <- c("MTB21_Lrg_BCP", "MTB21_Sml_BCP", "MTB21_MPFcom_BCP", "MTB09_Lrg_BCP", "MTB09_Sml_BCP", "MTB09_MPFcom_BCP")
  
  phyloseq.obj_norare.BCF.BCE.BCP <- merge_phyloseq(phyloseq.obj_norare.BCF.BCE, BCP)
  
  phyloseq.obj_norare.BCF.BCE.BCP 
  
}

############# FUNCTION TO MERGE SAMPLES AND CREATE THE BIOINFORMATICALLY COMBINED SAMPLES ######### 

Edesign.LM.pipeline <- function(phyloseq.obj, contrast.mat, design, R.size){
    
  # Transform non-orthogonal Contrast matrix
  newcontmatrix <- cbind(1, contrast.mat)
  # Create the the inverse of the transpose matrix  
  newcontmatrix <- solve(t(newcontmatrix))
  
  seeds <- c(NA,1000,2000,3000)
  out <- list()
  for(i in 1:4){
    #### Rarefying at even depth all the assigned and unassigned reads
    if(i == 1){
      #Remove the two samples with low read abundance
      sample_data(phyloseq.obj)@.Data[[1]] <- rownames(sample_data(phyloseq.obj))
      phyloseq.obj_allr_tmp <- subset_samples(phyloseq.obj, Sample_ID !="MTB11-B-01-PP" & Sample_ID !="MTB20-A-02-PP")
    }else{
      set.seed(seeds[i])  
      phyloseq.obj_allr_tmp <- rarefy_even_depth(phyloseq.obj, sample.size = R.size)
    }
    # Transform to PRESENCE/ABSENCE matrix 
    phyloseq.obj_allr_tmp_PA <- transform_sample_counts(phyloseq.obj_allr_tmp, function(abund) 1*(abund>0))
    
    OTU_table <- otu_table(phyloseq.obj_allr_tmp_PA)
    #### Analysis with OBSERVED RICHNESS 
    
    # Estimate observed richness
    phyloseq.obj_allr_tmp_PA_richness <- estimate_richness(OTU_table, measures ="Observed")
    phyloseq.obj_allr_tmp_PA_richness <-cbind(phyloseq.obj_allr_tmp_PA_richness, design)
    
    # Define the new contrast with the one associate in our data and run the linear models 
    contrasts(phyloseq.obj_allr_tmp_PA_richness$SFEP) <- newcontmatrix[,-1]
    
    #Run mixed model
    m <- lmer(Observed ~ SFEP+ (1|Smp), data = phyloseq.obj_allr_tmp_PA_richness)
    m.summary <- summary(m)  
    
    #### Permanova - Adonis with Jaccard dissim 
    PA <- t(OTU_table)
    rownames(PA) <- design$Sample_ID
    Jac <- vegdist(PA, method = "jaccard")
    
    treat.sum <- data.frame(model.matrix( ~ phyloseq.obj_allr_tmp_PA_richness$SFEP + phyloseq.obj_allr_tmp_PA_richness$Smp)[,-1])
    perm <- how(nperm = 999)
    setBlocks(perm) <- with(treat.sum, phyloseq.obj_allr_tmp_PA_richness$Smp)
    
    ad.mod <- adonis2(Jac ~ ., data = treat.sum, permutations = perm, by = "margin", contr.unordered = "contr.sum")
    
    #### Permadisp analysis 
    
    index <- apply(contrast.mat[,1:11], 2, function(x) x != 0)
    
    des.sub <- list()
    for(j in 1:ncol(index)){
      des <- droplevels(design[design$SFEP %in% levels(design$SFEP)[index[,j]],  c("Sample_ID","SFEP")])
      lvl <- data.frame(c.mat = contrast.mat[, j][index[,j]], ID = levels(design$SFEP)[index[,j]])
      des$lvl <- factor(lvl$c.mat[match(des$SFEP, lvl$ID)])
      des.sub[[j]] <- des
    }
    
    permdisp.l <- list()
    for(k in 1:length(des.sub)){
      tmp.des <- des.sub[[k]]
      tmp.mat <- as.matrix(Jac)[attr(Jac, "Labels") %in% as.character(tmp.des$Sample_ID), attr(Jac, "Labels") %in% as.character(tmp.des$Sample_ID)]
      b.disp <- betadisper(as.dist(tmp.mat), tmp.des$lvl, bias.adjust = TRUE)
      permdisp.l[[k]] <- b.disp
    }
    
    #permutest(betadisper(Jac, design$SFEP), permutations = perm, pairwise = TRUE)
    
    permdisp.anova <- lapply(permdisp.l, function(x) anova(x, permutations = perm))
    names(permdisp.anova) <- colnames(contrast.mat[,1:11])
    
    out[[i]] <- list(LM_richness = m.summary,  PERMANOVA_composition = ad.mod, PERMDISP_composition = permdisp.anova)
  }
  names(out) <- c("R0", paste("R", R.size,"_1", sep=""), paste("R", R.size,"_2", sep=""), paste("R", R.size,"_3", sep=""))
  out
}

#### Results NO Filtering ###### 

# Add the Sample_data to merge the corresponding samples
sample_data_tomerge <- read.csv(file = "Resources/sample_data_tomerge_edit.csv")

# Contrast matrix for the comparisons 
square_cm3 <- read.csv("Resources/square_cm3.csv")

# All samples used for the comparisons with its corresponding level  
Design_long_format_edit <- read.csv("Resources/Design_long_format_edit.csv")
Design_long_format_edit <- droplevels(Design_long_format_edit)
Design_long_format_edit$Smp <- as.factor(Design_long_format_edit$Smp)
Design_long_format_edit$SFEP <- as.factor(Design_long_format_edit$SFEP)
Design_long_format_edit

#### NOW RUN ALL ANAYSIS AGAIN TO MAKE SURE IT'S IN ORDER #### 
#### Results NO filtering ######

ASV_nofilter
ASV_nf.out  <- Edesign.merging(phyloseq.obj.raw = ASV_nofilter, sample.data = sample_data_tomerge)
ASV_nf.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_nf.out, contrast.mat = square_cm3, design = Design_long_format_edit, R.size = 4000)

ASV_nf_ass = subset_taxa(ASV_nofilter, Domain=="d__Eukaryota")
ASV_nf_ass.out  <- Edesign.merging(phyloseq.obj.raw = ASV_nf_ass, sample.data = sample_data_tomerge)
ASV_nf.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_nf_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

OTU_nofilter.out  <- Edesign.merging(phyloseq.obj.raw = OTU_nofilter, sample.data = sample_data_tomerge)
OTU_nf.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_nofilter.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

OTU_nf_ass = subset_taxa(OTU_nofilter, Domain=="d__Eukaryota")
OTU_nf_ass.out  <- Edesign.merging(phyloseq.obj.raw = OTU_nf_ass, sample.data = sample_data_tomerge)
OTU_nf.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_nf_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

#### Results LULU ###### 

ASV_lulu
ASV_lulu.out  <- Edesign.merging(phyloseq.obj.raw = ASV_lulu, sample.data = sample_data_tomerge)
ASV_lulu.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_lulu.out, contrast.mat = square_cm3, design = Design_long_format_edit, R.size = 4000)

ASV_lulu_ass = subset_taxa(ASV_lulu, Domain=="d__Eukaryota")
ASV_lulu_ass.out  <- Edesign.merging(phyloseq.obj.raw = ASV_lulu_ass, sample.data = sample_data_tomerge)
ASV_lulu.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_lulu_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

OTU_lulu
OTU_lulu.out <- Edesign.merging(phyloseq.obj.raw = OTU_lulu, sample.data = sample_data_tomerge)
OTU_lulu.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_lulu.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

OTU_lulu_ass = subset_taxa(OTU_lulu, Domain=="d__Eukaryota")
OTU_lulu.ass.out <- Edesign.merging(phyloseq.obj.raw = OTU_lulu_ass, sample.data = sample_data_tomerge)
OTU_lulu.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_lulu.ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

#### Results Filtering 0.003% ###### 

ASV_f1
ASV_f1.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f1, sample.data = sample_data_tomerge)
ASV_f1.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f1.out, contrast.mat = square_cm3, design = Design_long_format_edit, R.size = 4000)

ASV_f1_ass = subset_taxa(ASV_f1, Domain=="d__Eukaryota")
ASV_f1_ass.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f1_ass, sample.data = sample_data_tomerge)
ASV_f1.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f1_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

OTU_f1
OTU_f1.out <- Edesign.merging(phyloseq.obj.raw = OTU_f1, sample.data = sample_data_tomerge)
OTU_f1.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f1.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

OTU_f1_ass = subset_taxa(OTU_f1, Domain=="d__Eukaryota")
OTU_f1.ass.out <- Edesign.merging(phyloseq.obj.raw = OTU_f1_ass, sample.data = sample_data_tomerge)
OTU_f1.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f1.ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

#### Results Filtering 0.01% ###### 

ASV_f2
ASV_f2.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f2, sample.data = sample_data_tomerge)
ASV_f2.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f2.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

ASV_f2_ass = subset_taxa(ASV_f2, Domain=="d__Eukaryota")
ASV_f2_ass.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f2_ass, sample.data = sample_data_tomerge)
ASV_f2.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f2_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 2000)

OTU_f2
OTU_f2.out <- Edesign.merging(phyloseq.obj.raw = OTU_f2, sample.data = sample_data_tomerge)
OTU_f2.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f2.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

OTU_f2_ass = subset_taxa(OTU_f2, Domain=="d__Eukaryota")
OTU_f2_ass.out <- Edesign.merging(phyloseq.obj.raw = OTU_f2_ass, sample.data = sample_data_tomerge)
OTU_f2.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f2_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 1500)

#### Results Filtering 0.05% ###### 

ASV_f3
ASV_f3.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f3, sample.data = sample_data_tomerge)
ASV_f3.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = ASV_f3.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

ASV_f3_ass = subset_taxa(ASV_f3, Domain=="d__Eukaryota")
ASV_f3_ass.out  <- Edesign.merging(phyloseq.obj.raw = ASV_f3_ass, sample.data = sample_data_tomerge)
ASV_f3.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f3_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 1500)

OTU_f3
OTU_f3.out  <- Edesign.merging(phyloseq.obj.raw = OTU_f3, sample.data = sample_data_tomerge)
OTU_f3.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f3.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 4000)

OTU_f3_ass = subset_taxa(OTU_f3, Domain=="d__Eukaryota")
OTU_f3_ass.out <- Edesign.merging(phyloseq.obj.raw = OTU_f3_ass, sample.data = sample_data_tomerge)
OTU_f3.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU_f3_ass.out, contrast.mat = square_cm2, design = Design_longformat_merge, R.size = 1500)

save(ASV,ASVtable_tsc,OTUtable_tsc,tax_table_newtaxonomy_8ranks, reference_seqs0, sampledata, matchlistASV, matchlistOTU, ASV_nofilter, OTU_nofilter, ASV_lulu, OTU_lulu, ASV_f1, OTU_f1, ASV_f2, OTU_f2, ASV_f3, OTU_f3, sample_data_tomerge, square_cm3, Design_long_format_edit, file = "Edesign_github.Rdata")

