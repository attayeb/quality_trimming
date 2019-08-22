require(phyloseq)
library(stringr)
library(ade4)
library(ggplot2)
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)
library(vegan)
library(magrittr)

distance <- function(x, y){
    return(sqrt(sum(unlist(lapply(x-y, function(z) z^2)))))
}

phy = import_biom('biom_table.biom', treefilename = "rep_set.tre")
samples = c("-01", "-02", "-03", "-04", "-05", "-06", "-07", "-08", "-09", "-10")
pco_sample <- function(sample, phylo, plotting=FALSE, method='bray'){
    sphylo <- prune_samples(str_detect(sample_names(phylo), sample), phylo)
    flist = filterfun(kOverA(5, 10))
    filtered_otu_table <- otu_table(filter_taxa(s, flist, prune = TRUE))
    golden_standard <- apply(filtered_otu_table[,str_detect(colnames(filtered_otu_table) , "perfect")], 1 , mean)
    golden_standard <- data.frame(Var2 = "std", value=golden_standard)
    
    golden_standard$Var1 <- rownames(golden_standard)
    rownames(golden_standard) <- NULL
    filtered_otu_table <- filtered_otu_table[,!str_detect(colnames(filtered_otu_table) , "perfect")]
    
    melted_otu_table <- reshape2::melt(filtered_otu_table)
    #print(golden_standard)
    melted_otu_table <- rbind(melted_otu_table, golden_standard)
    
    casted_otu_table <- casted_otu_tableast(data=melted_otu_table, Var1~Var2)
    
    rownames(casted_otu_table) <- casted_otu_table$Var1
    casted_otu_table$Var1 <- NULL
    
    if (method== "unifrac") { distance_matrix <- UniFrac(t(casted_otu_table)) } else {
    distance_matrix <- vegdist(t(casted_otu_table), method)
        }
    distance_matrix[is.na(distance_matrix)] <- 1
    pco <- NULL
    pco <- dudi.pco(distance_matrix, nf=12, scannf = FALSE)
    return(pco)
}

pco_list <- list()
for (sample in samples){
    print(sample)
    pcoa <- plot_pco_sample(sample, phy12, method = "jaccard")
    euc_distance <- as.data.frame(apply(pcoa$li, 1, distance, pcoa$li[c('std'),]))
    rownames(euc_distance) -> euc_distance$sample
    colnames(euc_distance)[1] <- "distance"
    pco_list[[sample]] <- euc_distance
}