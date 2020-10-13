library(phyloseq)
library(decontam)
library(lme4)
library(lmerTest)
library(vegan)

#################################################################################
#### This function start with the 96 samples (including negative controls and other samples that are not part of this Experimental design, extract contaminants, merge the appropriate samples to create the Bioinformatically combined (BC) samples and subset only the samples used in the downstream analysis (linear models))

Edesign.cleaning <- function(phyloseq.obj.raw, sample.data){

# Identify and extract contaminants. The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants ###

sample_data(phyloseq.obj.raw)$is.neg <- sample_data(phyloseq.obj.raw)$Sample_or_Control == "Control"
contam.comb.all <- isContaminant(phyloseq.obj.raw, method="combined", neg="is.neg", threshold=0.5, conc="Concentration"  )
phyloseq.obj.nocon <- prune_taxa(!contam.comb.all$contaminant, phyloseq.obj.raw)

#### Merge samples  ###########

########## Create ASV table AND Subset ASV_table with the Experimental design samples + Controls 
phyloseq.obj <- subset_samples(phyloseq.obj.nocon, Edesign=="yes")
phyloseq.obj.nc <- subset_samples(phyloseq.obj, Sample_or_Control=="True Sample")

########## Merging samples that will be compared in our experimental design

#### Merge samples first and then rarefy!!!!
# Create a phyloseq object from the rarefy otu_table, sample data that allow us to merge by factor (sample_data_wo_raresmp)
phyloseq.obj.nc.norare <- phyloseq(otu_table(phyloseq.obj.nc), sample_data(data.frame(sample.data, row.names = sample_names(phyloseq.obj.nc))), tax_table(phyloseq.obj.nc), refseq(phyloseq.obj.nc)) #CHECK THIS STEP AFTER EDITING PHYLOSEQ OBJECT

## A+B (subset mA.B, merge samples by Sample and merge with the previous phyloseq) 
A.B <- subset_samples(phyloseq.obj.nc.norare, mAB=="mAB")
mA.B <- merge_samples(A.B,group = "Smp", fun = sum)
sample_names(mA.B) <- c("MTB09_mA.B", "MTB11_mA.B", "MTB21_mA.B")

phyloseq.obj_norare.AB <- merge_phyloseq(phyloseq.obj.nc.norare, mA.B)

## E1+E2 (subset mE1.E2, merge samples by Sample and merge with the previous phyloseq) 
E1.E2 <- subset_samples(phyloseq.obj.nc.norare, mE1.E2=="mE1.E2")
mE1.E2 <- merge_samples(E1.E2,group = "Smp", fun = sum)
sample_names(mE1.E2) <- c("MTB05_mE1.E2", "MTB09_mE1.E2", "MTB11_mE1.E2", "MTB14_mE1.E2", "MTB20_mE1.E2", "MTB21_mE1.E2", "MTB27_mE1.E2", "MTB30_mE1.E2")

phyloseq.obj_norare.AB.E12 <-merge_phyloseq(phyloseq.obj_norare.AB, mE1.E2)

## P1+P2 (subset mE1.E2, merge samples by Sample and merge with the previous phyloseq)
P1.P2 <- subset_samples(phyloseq.obj.nc.norare, mP1.P2=="mP1.P2")
mP1.P2 <- merge_samples(P1.P2,group = "SF", fun = sum)
sample_names(mP1.P2) <- c("MTB21_A_mP1.P2", "MTB21_B_mP1.P2", "MTB21_NS_mP1.P2", "MTB09_A_mP1.P2", "MTB09_B_mP1.P2", "MTB09_NS_mP1.P2")

phyloseq.obj_norare.AB.E12.P12 <- merge_phyloseq(phyloseq.obj_norare.AB.E12, mP1.P2)

phyloseq.obj_norare.AB.E12.P12

}

#################################################################################
##### Edesign.LM.pipeline function uses the contrast matrix we used for our tests and run the linear models with the rarefied presence/absence matrix 

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
    #### Analysis with OBSERVED RICHNESS ####
    
    # Estimate observed richness
    phyloseq.obj_allr_tmp_PA_richness <- estimate_richness(OTU_table, measures ="Observed")
    phyloseq.obj_allr_tmp_PA_richness <-cbind(phyloseq.obj_allr_tmp_PA_richness, design)
    
    # Define the new contrast with the one associate in our data and run the linear models 
    contrasts(phyloseq.obj_allr_tmp_PA_richness$SFEP) <- newcontmatrix[,-1]
    
    #Run mixed model
    m <- lmer(Observed ~ SFEP+ (1|Smp), data = phyloseq.obj_allr_tmp_PA_richness)
    m.summary <- summary(m)  
    
    #### Analysis with Jaccard dissim ####
    PA <- t(OTU_table)
    
    Jac <- vegdist(PA, method = "jaccard")
    
    treat.sum <- data.frame(model.matrix( ~ phyloseq.obj_allr_tmp_PA_richness$SFEP + phyloseq.obj_allr_tmp_PA_richness$Smp)[,-1])
    perm <- how(nperm = 999)
    setBlocks(perm) <- with(treat.sum, phyloseq.obj_allr_tmp_PA_richness$Smp)
    
    ad.mod <- adonis2(Jac ~ ., data = treat.sum[,-25], permutations = perm, by = "margin", contr.unordered = "contr.sum")
    out[[i]]<- list(LM_richness = m.summary,  PERMANOVA_composition = ad.mod)
  }
  names(out) <- c("R0", paste("R", R.size,"_1", sep=""), paste("R", R.size,"_2", sep=""), paste("R", R.size,"_3", sep=""))
  out
}



###### All the R.data that we need to get the Results for the 16 tables ((x4) 3 rarefactions replicates and no-rarefied))

load("~/Desktop/Edesign manuscript/Script and R data/Edesign/PS_objects_Edesign.Rdata")
load("~/Desktop/Edesign manuscript/Script and R data/Edesign/Sample_data_tomerge63.Rdata")

#### Results NO Filtering ###### 

ASV_nf.out  <- Edesign.cleaning(phyloseq.obj.raw = ASV_nf, sample.data = sample_data_tomerge)

ASV_nf.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_nf.out, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

ASV_nf.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 2500)

OTU_nf.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU00_allr, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

OTU_nf.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU00_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 2500)


#### Results Filtering 0.003% ###### 

ASV_f39.out  <- Edesign.cleaning(phyloseq.obj.raw = ASV_f39, sample.data = sample_data_tomerge)

ASV_f39.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f39.out, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

ASV_f39.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV39_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 2500)

OTU_f39.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU39, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

OTU_f39.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU39_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 2500)

#### Results Filtering 0.01% ###### 

ASV_f130.out  <- Edesign.cleaning(phyloseq.obj.raw = ASV_f130, sample.data = sample_data_tomerge)

ASV_f130.LM.allr.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV_f130.out, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

ASV_f130.LM.ass.out  <- Edesign.LM.pipeline(phyloseq.obj = ASV130_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 2000)

OTU_f130.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU130, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

OTU_f130.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU130_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 1500)

#### Results Filtering 0.05% ###### 

ASV_f652.out  <- Edesign.cleaning(phyloseq.obj.raw = ASV_f652, sample.data = sample_data_tomerge652)

ASV_f652.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = ASV_f652.out, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

ASV_f652.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = ASV652_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 1800)

OTU_f652.LM.allr.out <- Edesign.LM.pipeline(phyloseq.obj = OTU652, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 4000)

OTU_f652.LM.ass.out <- Edesign.LM.pipeline(phyloseq.obj = OTU652_Ed_ass, contrast.mat = square_cm, design = Design_longformat_merge, R.size = 1500)
