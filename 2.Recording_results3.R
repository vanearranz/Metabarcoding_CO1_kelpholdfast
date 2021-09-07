########################################################
## Script to record the results from Edesign_function_results.R. Figure in Supporting material A2 and Summary results in Figure 2 in the manuscrript
########################################################

library(ggplot2)
library(ggh4x)

load("Resources/Edesign_results.Rdata")

results.list <- list(ASV_nf.LM.allr.out, ASV_lulu.LM.allr.out, ASV_f1.LM.allr.out, ASV_f2.LM.allr.out, ASV_f3.LM.allr.out, ASV_nf.LM.ass.out, ASV_lulu.LM.ass.out, ASV_f1.LM.ass.out, ASV_f2.LM.ass.out, ASV_f3.LM.ass.out, OTU_nf.LM.allr.out, OTU_lulu.LM.allr.out, OTU_f1.LM.allr.out, OTU_f2.LM.allr.out, OTU_f3.LM.allr.out, OTU_nf.LM.ass.out, OTU_lulu.LM.ass.out, OTU_f1.LM.ass.out, OTU_f2.LM.ass.out, OTU_f3.LM.ass.out)

names(results.list) <- paste(rep(c("ASV", "OTU"), each = 10), rep(c("All", "Taxonomically assigned"), each = 5, times = 2), rep(c("NO", "LULU", "0.003%", "0.01%", "0.05%"), times = 4), sep= "_")

tile.sum.all <- list()
for(j in 1:length(results.list)){
  tmp <- results.list[[j]]
  tile.sum.tmp <- list()
    for(i in 1:length(tmp)){
  
    lm.col <- ifelse(round(tmp[[i]]$LM_richness$coefficients[2:12,5],3) <0.001, "A", ifelse(round(tmp[[i]]$LM_richness$coefficients[2:12,5],3) <0.01, "B", ifelse(round(tmp[[i]]$LM_richness$coefficients[2:12,5],3) <0.05, "C", "D")))
    
    adonis.col <- ifelse(tmp[[i]]$PERMANOVA_composition$`Pr(>F)`[1:11] <0.001, "A", ifelse(tmp[[i]]$PERMANOVA_composition$`Pr(>F)`[1:11] <0.01, "B", ifelse(tmp[[i]]$PERMANOVA_composition$`Pr(>F)`[1:11] <0.05, "C", "D")))
  
    permdisp.p <- unlist(lapply(tmp[[i]]$PERMDISP_composition, function(x) x$`Pr(>F)`[1]))
  
    permdisp.col <- ifelse(permdisp.p <0.001, "A", ifelse(permdisp.p <0.01, "B", ifelse(permdisp.p <0.05, "C", "D")))
    
    tile.df <- data.frame(R.seed = c("R0","R1","R2","R3")[i],  Model = factor(rep(c("Richness", "Composition", "Dispersions"), each = 11), levels = (c("Richness", "Composition", "Dispersions"))), contrast = factor(rep(names(permdisp.col),3), levels = (names(permdisp.col))), Sig. = c(lm.col, adonis.col, permdisp.col))  
    tile.sum.tmp[[i]] <- tile.df
    }
  index <- unlist(strsplit(names(results.list)[j], "_"))
tile.sum.all[[j]] <- data.frame(Cluster = index[1], Assigned = index[2], Curation = index[3], do.call("rbind", tile.sum.tmp))
}

tile.sum.all.p <- do.call("rbind", tile.sum.all)
tile.sum.all.p$Cluster <- factor(tile.sum.all.p$Cluster, levels = c("ASV", "OTU"))
tile.sum.all.p$Assigned <- factor(tile.sum.all.p$Assigned, levels = c("All", "Taxonomically assigned"))
tile.sum.all.p$Curation <- factor(tile.sum.all.p$Curation, levels = c("NO", "LULU", "0.003%", "0.01%", "0.05%"))

SM_A2 <- ggplot(tile.sum.all.p)+
  geom_tile(aes(y = contrast, x = R.seed, fill = Sig.)) +
  scale_fill_manual(values = c("#204C7B", "#3278BE", "#C8E3F2", "white"), labels = c("<0.001", "<0.01", "<0.05", ">0.05")) +
  facet_nested(Model + contrast ~ Cluster + Assigned + Curation + R.seed, scales = "free") +
  scale_x_discrete(position = "top") +
  theme(strip.text.y = element_text(angle = 0, size = 6), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.y.left = element_blank(), axis.title.y.left = element_blank(), strip.text.x = element_text(face="bold", size=6), axis.title.x=element_blank(), axis.text.x=element_blank(), panel.spacing.y = unit(0.05, "lines")) 

#pdf("Figure A2.pdf", width = 16, height = 11, family = "Courier", useDingbats=FALSE)
SM_A2
#dev.off() 
