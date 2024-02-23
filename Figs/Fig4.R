library(ggplot2)
library(stringr)
library(reshape2)
library(tidyr)
library(gridExtra)
library(tidyverse)
library(pheatmap)
library(grid)
library(cowplot)
library(RColorBrewer)

counts<-read.csv("Rice_phosphosite_matrix_binomial_peptidoform_0.05_w_protein-pos_scores_FLR_SNP_reps_expand_per_PTM_oryzabase_RAP-DB_MSU_freq_var_all_no_pA.csv")

#a) bar chart counts SNP category
cat_counts<-subset(counts,select=c("PTM_site_category"))
cat_counts<-melt(table(cat_counts))
a<-ggplot(cat_counts,aes(y=value, x=cat_counts))+ 
  geom_bar(position='dodge', stat='identity')+
  ylab("Count of Sites")+
  xlab("Category")+
  ggtitle("a)")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )



# b) cat1 heatmap, sort by absolute diff between jap family vs indica family allel freq
#x=4 families
#y=cat1 sites with gene label (gene_site)
#heat=allele freq
counts_cat1<-subset(counts[counts$PTM_site_category == "1",])
counts_cat1<-counts_cat1[!is.na(counts_cat1$PTM_site_category),]
counts_cat1<-subset(counts_cat1[counts_cat1$Annotated. == "True",], 
                    select=c("jap_family_ref_freq","ind_family_ref_freq","admx_family_ref_freq","aus_family_ref_freq",
                             "Total_diff","Jap_Ind_Diff","PTM_pos","Oryzabase","RAP_genes","MSU_genes","PTM_site_category"))
counts_cat1$Gene_Pos<-paste(counts_cat1$Oryzabase, counts_cat1$RAP_genes, counts_cat1$MSU_genes, sep=":")
counts_cat1$Gene_Pos<-gsub("_;_","",counts_cat1$Gene_Pos)
counts_cat1$Gene_Pos<-gsub("-;-","",counts_cat1$Gene_Pos)
counts_cat1$Gene_Pos<-gsub(":;","",counts_cat1$Gene_Pos)
counts_cat1$Gene_Pos<-gsub("::","",counts_cat1$Gene_Pos)
counts_cat1$Gene_Pos<-paste(counts_cat1$Gene_Pos, counts_cat1$PTM_pos, sep="_")
counts_cat1<-counts_cat1[!duplicated(counts_cat1$Gene_Pos),]
rownames(counts_cat1) <- counts_cat1$Gene_Pos
counts_cat1<-counts_cat1[order(counts_cat1$Jap_Ind_Diff),]
#filter for total_diff > 0.25
counts_cat1<-subset(counts_cat1[counts_cat1$Total_diff>0.05,])
counts_cat1<-counts_cat1[, c("jap_family_ref_freq","ind_family_ref_freq","admx_family_ref_freq","aus_family_ref_freq")]

counts_cat1[is.na(counts_cat1)]<-0

counts_cat1 <- as.matrix(counts_cat1)

cluster_no <- 5
col <-colorRampPalette(c("white","navy"))(256)
heatmap <- pheatmap(counts_cat1, cutree_rows = cluster_no,
                 show_rownames=FALSE,
                 clustering_distance_rows="euclidean",
                 clustering_distance_cols="euclidean",
                 clustering_method="complete",
                 main="c)                                                                                                                                                   ",
                 color=col,
                 labels_col=c("Japonica","Indica", "Admx", "Aus"))
                 #,legend_breaks = c(0, 25, 50, 75,100),legend=TRUE,
                 #legend_labels = c("0%", "25%", "50%", "75%","100%"))
#b<-grobTree(heatmap)

cluster_names <- cutree(heatmap$tree_row, k = cluster_no)
cluster_names[1:44]<-""
cluster_names[47:86]<-""
cluster_names[88:94]<-""
cluster_names[96:108]<-""
cluster_names[110:125]<-""

cluster_names[45]<-"Cluster 4"
cluster_names[46]<-"Cluster 5"
cluster_names[87]<-"Cluster 1"
cluster_names[95]<-"Cluster 2"
cluster_names[109]<-"Cluster 3"

heatmap <- pheatmap(counts_cat1, cutree_rows = cluster_no,
                    show_rownames=TRUE,
                    clustering_distance_rows="euclidean",
                    clustering_distance_cols="euclidean",
                    clustering_method="complete",
                    main="c)                                                                                                                                                   ",
                    color=col,
                    labels_col=c("Japonica","Indica", "Admx", "Aus"),
                    labels_row = cluster_names)


#extract data into specified cluster (note: the assigned numbers of clusters may not match the order in heatmap)
clusters <- as.data.frame(cbind(counts_cat1, cluster = cutree(heatmap$tree_row, k = cluster_no)))
clusters$cluster[clusters$cluster == 1] <- 'Cluster 5'
clusters$cluster[clusters$cluster == 2] <- 'Cluster 4'
clusters$cluster[clusters$cluster == 3] <- 'Cluster 1'
clusters$cluster[clusters$cluster == 4] <- 'Cluster 2'
clusters$cluster[clusters$cluster == 5] <- 'Cluster 3'

#list all row names in the order of heatmap
all_rows <- rownames(counts_cat1[heatmap$tree_row[["order"]],])

#c) cat1, counts aa mutated from ST vs normalised aa count
counts_cat1<-subset(counts[counts$PTM_site_category == "1",])
counts_cat1<-counts_cat1[!is.na(counts_cat1$PTM_site_category),]
counts_cat1[c("SNP_before", "SNP_after")]<-str_split_fixed(counts_cat1$SNP_res,"->",2)
SNP_after_counts<-subset(counts_cat1,select=c("SNP_after"))
## keep * values
SNP_after_counts$SNP_after<- stringr::str_replace(SNP_after_counts$SNP_after, '\\*', '.')
SNP_after_counts<-separate_rows(SNP_after_counts,SNP_after,convert=TRUE)
SNP_after_counts$SNP_after<- stringr::str_replace(SNP_after_counts$SNP_after, '\\.', '*')
SNP_after_counts<-melt(table(SNP_after_counts))
count_sites <- sum(SNP_after_counts$value)

aa_count<-subset(counts,select=c("Peptide"))
aa_count$Peptide<-gsub("\\-.*","",aa_count$Peptide)
aa_count$Peptide<-gsub("(?<=.)(?=.)", ",", aa_count$Peptide, perl = TRUE)
aa_count<-separate_rows(aa_count,Peptide,convert=TRUE)
aa_count<-melt(table(aa_count))
aa_count_sites<- sum(aa_count$value)
aa_count$value<-(aa_count$value/aa_count_sites)*count_sites

colnames(SNP_after_counts)[1]<-"residue"
SNP_after_counts$Distribution <- "SNP mutated residues"
colnames(aa_count)[1]<-"residue"
aa_count$Distribution <- "background amino acid distribution"

all_counts <- rbind(SNP_after_counts,aa_count)
c<-ggplot(all_counts, aes(residue, value, fill=Distribution))+
  geom_bar(stat="identity", position=position_dodge(preserve="single"))+
  ggtitle("b)")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

lm <- rbind(c(1,2),
            c(3,3),
            c(3,3))
p4<-grid.arrange(grobs=list(a,c,heatmap[[4]]), layout_matrix=lm)

