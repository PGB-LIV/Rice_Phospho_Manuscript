require(rmotifx)
require(clusterProfiler)
library(stringr)
require(ggplot2)
require(ggseqlogo)
require(gridExtra)
require(Rmisc)
library(dplyr)
library(tidyr)
library(pRoloc)
library(cowplot)

wd="All_datasets"
setwd(wd)
dir.create(file.path(wd, "Motif_analysis"))

# Read in sequences
background<-read.csv(file="motif_background.txt", header=FALSE)
bg.seqs = as.character(background$V1)
bg_list<-unique(as.character(background$V3))

source('mapping_and_enrichment-main/mapping_and_enrichment-main/uniprot_dat_reading.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/mapping_with_db_ortho.R')
##Make sure right species pointing here!
source('mapping_and_enrichment-main/mapping_and_enrichment-main/mapping_with_files.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/uniprot_selected_term_types.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/utils.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/enrichment.R')


#ClusterProfiler background - all phosphoproteins
cluster_background<-read.csv(file="All_phospoproteins_background.txt", header=FALSE)
cluster_bg_list<-unique(as.character(cluster_background$V1))


file_list<-c("gsb_motif_seqs.txt","gs_motif_seqs.txt","gold_motif_seqs.txt")
for (f in file_list){
  foreground<-read.csv(file=f,header=FALSE)
  fg.seqs = as.character(foreground$V1)
  print(f)
  #find enriched motifs
  if (f=="S_motif_seqs.txt"){
    mot = motifx(fg.seqs, bg.seqs, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
  } else if (f=="T_motif_seqs.txt"){
    mot = motifx(fg.seqs, bg.seqs, central.res = 'T', min.seqs = 20, pval.cutoff = 1e-6)
  } else if (f=="Y_motif_seqs.txt"){
    mot = motifx(fg.seqs, bg.seqs, central.res = 'Y', min.seqs = 20, pval.cutoff = 1e-6)
  } else{
    mot = motifx(fg.seqs, bg.seqs, central.res = 'ST', min.seqs = 20)
  }
  
  write.csv(mot, file = paste0("Motif_analysis/All_motifs_",f,".csv"))
  
  #create GO dataframe
  GO<-data.frame()
  
  #motif seqlogo of enriched motif sequences
  pos<-c("-7","-6","-5","-4","-3","-2","-1","p","1","2","3","4","5","6","7")
  break_list<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  plots<-list()
  plots2<-list()
  counter=1
  mid=round(length(mot$motif)/2)
  for (val in mot$motif){
    m<-grep(val,fg.seqs,value=TRUE)
    motif_plot=ggplot()+geom_logo(m, method="probability")+scale_x_continuous(labels=pos,breaks=break_list)+ggtitle(val)+theme(legend.position = "none")
    plots[[counter]]<-motif_plot
    
    #Pathway enrichment analysis - ClusterProfiler
    
    # All phosphosites
    # Define foreground and background gene lists.
    # The foreground list should be contained within the background list.
    foreground_motif<-foreground[grep(val,foreground$V1),]
    fg_list<-unique(as.character(foreground_motif$V3))
    
    background_terms <- map_using_uniprot_files(input_id = 'UniProtKB-AC', input_values = cluster_bg_list, output_id = 'GO')
    background_terms <- background_terms[,c(2,1)]
    background_terms$names <- goIdToTerm(background_terms$GO, names = TRUE, keepNA = TRUE)
    
    enrichment <- enricher(
      fg_list,
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      universe=as.character(cluster_bg_list),
      minGSSize = 3,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = background_terms[c(1,2)],
      TERM2NAME = background_terms[c(1,3)]
    )
    enrichment
    if (!is.null(enrichment)){
      if (enrichment@result[["p.adjust"]][1]<enrichment@pvalueCutoff){
        GO_temp<-as.data.frame(enrichment@result)
        GO_temp$motif<-val
        GO<-rbind(GO,GO_temp)
        ep<-clusterProfiler::dotplot(enrichment, title=val)
        try(plots2[[counter]]<-ep)
      }else{
        ep<-ggplot() + ggtitle(val)
        try(plots2[[counter]]<-ep)
      }
      
    }else{
      ep<-ggplot() + ggtitle(val)
      try(plots2[[counter]]<-ep)
    }
    counter = counter+1
  
  write.csv(GO, file = paste0("Motif_analysis/CP_",f,".csv"))

  p2<-grid.arrange(grobs=plots)
  ggsave(file=paste0("Motif_analysis/Motif_all_",f,".png"),
         p2, dpi=300, width=8000, height=6500, units="px")

  p2<-grid.arrange(grobs=plots2)
  ggsave(file=paste0("Motif_analysis/CP_all_",f,".png"),
         p2, dpi=300, width=10000, height=12000, units="px")
  }
}



