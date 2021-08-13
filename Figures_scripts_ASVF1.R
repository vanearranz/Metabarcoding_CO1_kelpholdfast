########################################################
### Script to Plot Figures 3 and 4 in manuscript "Metabarcoding hyperdiverse marine communities in temperate kelp forests: an experimental approach to inform future studies." 
########################################################

library(phyloseq)
library(lme4)
library(lmerTest)
library(vegan)
library(emmeans)
library(ggplot2)
library(ggthemes)

load("~/Resources/Edesign.Rdata")
load("~/Resources/Edesign_results.Rdata")

### Figure 3 ####
## Considering ASV lulu CONSIDERING ALL READS 

ASV_f1.out

sample_data(ASV_f1.out)@.Data[[1]] <- rownames(sample_data(ASV_f1.out))
ASV_f1.out_tmp <- subset_samples(ASV_f1.out, Sample_ID!="MTB11-B-01-PP" & Sample_ID !="MTB20-A-02-PP")
  
#### R4000 : Rarefying at even depth all the assigned and unassigned reads 
set.seed(1000)
ASV_f1.out_tmp_R4000 = rarefy_even_depth(ASV_f1.out_tmp, sample.size = 4000)

# Transform to PRESENCE/ABSENCE matrix 
ASV_f1.out_tmp_R4000_PA <- transform_sample_counts(ASV_f1.out_tmp_R4000, function(abund) 1*(abund>0))

sample_data(ASV_f1.out_tmp)$SFEP <- Design_long_format_edit$SFEP

# For non-orthogonal contrast we have to transpose the inverse (I already have my square matrix because I transformed it previously))
newcontmatrix <- cbind(1, square_cm4)
# Create the the inverse of the transpose matrix  
newcontmatrix <- solve(t(newcontmatrix))
newcontmatrix
round(newcontmatrix, 2)

# Estimate observed richness
ASV_f1.out_tmp_R4000_PA_richness <- estimate_richness(otu_table(ASV_f1.out_tmp_R4000_PA), measures ="Observed")
ASV_f1.out_tmp_R4000_PA_richness<-cbind(ASV_f1.out_tmp_R4000_PA_richness, Design_long_format_edit)

#### Run the lm analysis for each of the richness estimates 

## Define where the diversity estimate is 
y <- ASV_f1.out_tmp_R4000_PA_richness$Observed

# Define the new contrast with the one associate in our data and run the linear models 

# Transform non-orthogonal Contrast matrix
newcontmatrix <- cbind(1, square_cm4)
# Create the the inverse of the transpose matrix  
newcontmatrix <- solve(t(newcontmatrix))


contrasts(ASV_f1.out_tmp_R4000_PA_richness$SFEP) <- newcontmatrix[,-1]
m <- lmer(y ~ SFEP + (1|Smp), data = ASV_f1.out_tmp_R4000_PA_richness)
summary(m)

#Get the group means to plot 

emm <- emmeans(m, "SFEP")

###Get the contrast that we had previously done
design <- square_cm4[,1:11]

contrast.names <- unlist(sapply(colnames(design), function(x) strsplit(x, "vs")))

LC.list <- list()
for(i in 1:dim(design)[2]){
  tmp <- design[,i]
  index <- ifelse(tmp == 0, -999, tmp)
  mm <- model.matrix(~ as.factor(index))[,-1]
  colnames(mm) <- rev(unlist(strsplit(colnames(design)[i], "vs")))
  LC.list[[i]] <- mm
}
LC.df <- do.call("cbind", LC.list)
c.sum <- colSums(LC.df)

LC.prop <- LC.df
for(i in 1:dim(LC.df)[2]){
  LC.prop[,i] <- LC.df[,i]/c.sum[i]
}

# Means for each factor of the contrasts seed.1
LF <- contrast(emm, list(LC.prop))
confint(LF)

LF.df <- as.data.frame(LF)

#We plot only one seed LF.df
#ggplot will plot strong variables as factors and factors are alphabetically ordered these lines tell it to leave the order as its is
SP.o <- c("Sml", "Lrg", "MPFpcr", "BAF","BCF", "MPFext","MPFcom", "EA", "EB", "MPEext", "BAE", "BCE", "PX", "PY", "MPPpcr", "BAP", "BCP")

LF.df$contrast <- as.character(LF.df$contrast)
LF.df$contrast <- factor(LF.df$contrast, levels=SP.o)

#To set up the panels
LF.df$panel <- factor(rep(c("Sample Preparation", "DNA Extraction", "PCR Amplification"), times = c(10,6,6)), levels = c("Sample Preparation", "DNA Extraction", "PCR Amplification"))

plot.means <- ggplot(LF.df) + theme(panel.border = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(aes(x = contrast, y = estimate))  +
  geom_linerange(aes(x = contrast, ymin = estimate - SE, ymax = estimate + SE)) +
  facet_grid(cols = vars(panel), scales = "free_x") +
  labs(x="Experimental Treatment", y = "Observed richness")
  
plot.means  

### Figure 4 ####

#### Spider plot for the Small versus Large comparison 

ASV_f1.out_tmp_R4000_PA

ASV_f1.out_tmp_R4000_PA_table <- as.data.frame(otu_table(ASV_f1.out_tmp_R4000_PA))
PA <- t(ASV_f1.out_tmp_R4000_PA_table)
rownames(PA) <- Design_long_format_edit$Sample_ID

Jac <- vegdist(PA, method = "jaccard")

treat.sum <- data.frame(model.matrix( ~ ASV_f1.out_tmp_R4000_PA_richness$SFEP + ASV_f1.out_tmp_R4000_PA_richness$Smp)[,-1])
AvsB <- ifelse(round(treat.sum[,1],1)<0, "Small", ifelse(round(treat.sum[,1],1)>0, "Large", NA))

beta.JacAvsB <- betadisper(Jac, AvsB)

beta.JacAvsB$group

coords <- data.frame(scores(beta.JacAvsB)$sites)
cen <- data.frame(scores(beta.JacAvsB)$centroids)

coord.hull.ls <- list()
for (i in levels(beta.JacAvsB$group)) {
  ch <- chull(coords[beta.JacAvsB$group == i, ])
  ch <- c(ch, ch[1])
  coord.hull.ls[[i]] <- data.frame(coords[beta.JacAvsB$group == i, ][ch,], Fraction = i)
}

smp.ls <- strsplit(rownames(coords), "-")
coords$Smp <- as.factor(unlist(lapply(smp.ls, function(x) paste(x[1]))))
coords$Fraction <- as.factor(unlist(lapply(smp.ls, function(x) x[2])))
coords$Fraction <- ifelse(coords$Fraction == "B", "Small", "Large")

coord.hull.df <- data.frame(do.call("rbind", coord.hull.ls))
eig.p <- (beta.JacAvsB$eig/sum(beta.JacAvsB$eig))*100

spider_AvsB <- ggplot() + theme_base()+
  geom_polygon(data = coord.hull.df, aes(x = PCoA1, y = PCoA2, group = Fraction, colour = NULL, fill = Fraction), alpha = 0.3) +
  geom_point(data = coords, aes(x = PCoA1, y = PCoA2, shape = Fraction, colour = Smp), size = 3) +
  geom_text(data = cen, aes(x = PCoA1, y = PCoA2, label = rownames(cen)))+
  scale_fill_manual(values = c("steelblue2", "firebrick2"), labels = c("Small", "Large"), guide = "none")+
  xlab(paste("PCoA1 ", "(",round(eig.p[1], digits = 1), "%)", sep = "")) + 
  ylab(paste("PCoA2 ", "(",round(eig.p[2], digits = 1), "%)", sep = "")) +
  theme(legend.position = "none")

spider_AvsB 

