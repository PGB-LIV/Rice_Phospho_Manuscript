require(rmotifx)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(cowplot)
require(ggplot2)
require(ggseqlogo)
require(gridExtra)
require(Rmisc)
library(grid)
require(clusterProfiler)
library(pRoloc)

wd="All_datasets"
setwd(wd)

source('mapping_and_enrichment-main/mapping_and_enrichment-main/uniprot_dat_reading.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/mapping_with_db_ortho.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/mapping_with_files.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/uniprot_selected_term_types.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/utils.R')
source('mapping_and_enrichment-main/mapping_and_enrichment-main/enrichment.R')

#ClusterProfiler background - all phosphoproteins
cluster_background<-read.csv(file="All_phospoproteins_background.txt", header=FALSE)
cluster_bg_list<-unique(as.character(cluster_background$V1))

######
protein_filter=100

# Read in sequences
background<-read.csv(file="motif_background.txt", header=FALSE)
bg.seqs = as.character(background$V1)
bg_list<-unique(as.character(background$V3))

file_list<-c("gsb_motif_seqs.txt")
for (f in file_list){
  foreground<-read.csv(file=f,header=FALSE)
  fg.seqs = as.character(foreground$V1)
  #find enriched motifs
  mot = motifx(fg.seqs, bg.seqs, central.res = 'ST', min.seqs = 20, pval.cutoff = 1e-6)
  
  unique_proteins<-list()
  
  for (val in mot$motif){
    
    if (file.exists(paste0("Motif_analysis_2/CP_",f,val,".csv"))){
  
      func_chart<-read.csv(paste0("Motif_analysis_2/CP_",f,val,".csv"),header=TRUE)
      if (nrow(func_chart)>0){
        
        #count no proteins associated with motif
        proteins<-func_chart %>% 
          separate_rows(geneID, sep="[//\r\n]+")
        unique_proteins<-c(unique_proteins,n_distinct(proteins$geneID))
        
        #filter func chart for count>3 and 1-FDR>0.6
        func_chart<-filter(func_chart, qvalue<=0.4)
        func_chart<-filter(func_chart, Count>3)
        
        #1-q_value
        temp<-data.frame(func_chart$Description,func_chart$qvalue)
        #FDR>1 - set to 1
        temp$func_chart.qvalue<-1-(temp$func_chart.qvalue)
        names(temp)[names(temp) == "func_chart.qvalue"] <- val
        names(temp)[names(temp) == "func_chart.Description"] <- "ID"
  
        if (val==mot$motif[1]){
          motifs<-temp
        }else{
          motifs<-dplyr::full_join(motifs, temp, by="ID")
        }
      }else{
        motifs[val]<-NA
        unique_proteins<-c(unique_proteins,"0")
        }
    }else{
      motifs[val]<-NA
      unique_proteins<-c(unique_proteins,"0")
    }
  }
  rownames(motifs) <- motifs$ID
  drops<-c('ID')
  motifs<-motifs[ , !(names(motifs) %in% drops)]
  motifs[is.na(motifs)] <- 0
  
  mot_col<-data.frame(as.numeric(matrix(unlist(unique_proteins))))
  colnames(mot_col)<-c("Unique Proteins")
  rownames(mot_col)<-mot$motif
  mot_col$`Site Count`<-mot$fg.matches
  
  #filter for unique proteins
  mot_col<-mot_col%>% filter(mot_col$`Unique Proteins`>protein_filter)
  filter<-as.list(rownames(mot_col))
  filter<-append(filter, "ID")
  motifs<-motifs[,(names(motifs) %in% filter)]
  
  #remove rows where all values are zero
  motifs<-motifs[rowSums(motifs[])>0,]
  
  motif_list<-colnames(motifs)
  mot<-filter(mot, mot$motif %in% motif_list)
  
  #motif seqlogo of enriched motif sequences
  pos<-c("-7","-6","-5","-4","-3","-2","-1","p","1","2","3","4","5","6","7")
  break_list<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  plots<-list()
  plots2<-list()
  counter=1
  counter2=1
  mid=round(length(mot$motif)/2)
  for (val in mot$motif){
    m<-grep(val,fg.seqs,value=TRUE)
    motif_plot=ggplot()+geom_logo(m, method="probability")+scale_x_continuous(labels=pos,breaks=break_list)+ggtitle(val)+theme(legend.position = "none")
    plots[[counter]]<-motif_plot
    counter=counter+1
    
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
      if (length(enrichment)>0){
          if (enrichment@result[["p.adjust"]][1]<enrichment@pvalueCutoff){
            ep<-clusterProfiler::dotplot(enrichment, title=val)
            try(plots2[[counter2]]<-ep)
            counter2 = counter2+1
          }
      }
    }
  }  
  p1<-grid.arrange(grobs=plots)
  
  p2<-grid.arrange(grobs=plots2)
  
  motifs<-as.matrix(motifs)
  
  cluster_no<-5
  heatmap <- pheatmap(motifs, 
                      show_rownames=FALSE,
                      cutree_rows = cluster_no,
                      clustering_distance_rows="euclidean",
                      clustering_distance_cols="euclidean",
                      clustering_method="complete",
                      annotation_col=mot_col,
                      color=rev(hcl.colors(4, "Heat2")),
                      background="transparent")
  
  #extract data into specified cluster (note: the assigned numbers of clusters may not match the order in heatmap)
  clusters <- cbind(motifs, cluster = cutree(heatmap$tree_row, k = cluster_no))
  
  #list all row names in the order of heatmap
  all_rows <- rownames(motifs[heatmap$tree_row[["order"]],])
}


#Annotations
cs <- scale_colour_gradientn(colours=rev(hcl.colors(4, "Heat2")), limits=c(0,1))

#Cluster 1 - [ST]P
func_list<-filter(as.data.frame(clusters), cluster==1)
func_list<-rownames(func_list)
ann1_temp<-filter(as.data.frame(clusters),cluster==1)
ann1<-data.frame(rownames(ann1_temp),ann1_temp$`.......[ST]P......`)
names(ann1)<-c("Description","1-FDR")
func_chart<-filter(func_chart, Description %in% func_list) 
ann1<-dplyr::full_join(ann1, subset(func_chart, select=c("Description", "GeneRatio","Count")), by="Description")
ann1$GeneRatio<-DOSE::parse_ratio(ann1$GeneRatio)
ann1<-filter(ann1, Count>3)
ann1<-filter(ann1, `1-FDR`>=0.6)
ggplot(ann1, aes(x=GeneRatio, y=Description, colour=`1-FDR`, size=Count))+geom_point()+
  cs+xlab("Enrichment")+ylab("")+ scale_size_continuous(limits=c(0,60),breaks=c(10,20,30,40,50))+theme_gray(base_size = 18)


#Cluster 4 - [ST]P.R
func_list<-filter(as.data.frame(clusters), cluster==4)
func_list<-rownames(func_list)
ann1_temp<-filter(as.data.frame(clusters),cluster==4)
ann1<-data.frame(rownames(ann1_temp),ann1_temp$`.......[ST]P.R....`)
names(ann1)<-c("Description","1-FDR")
func_chart<-read.csv(paste0("Motif_analysis_2/CP_gsb_motif_seqs.txt.......[ST]P.R.....csv"),header=TRUE)
func_chart<-filter(func_chart, Description %in% func_list) 
ann1<-dplyr::full_join(ann1, subset(func_chart, select=c("Description", "GeneRatio","Count")), by="Description")
ann1$GeneRatio<-DOSE::parse_ratio(ann1$GeneRatio)
ann1<-filter(ann1, Count>3)
ann1<-filter(ann1, `1-FDR`>=0.6)
ggplot(ann1, aes(x=GeneRatio, y=Description, colour=`1-FDR`, size=Count,na.rm = TRUE))+
  geom_point()+cs+xlab("Enrichment")+ylab("")+ 
  scale_size_continuous(limits=c(0,60),breaks=c(10,20,30,40,50))+
  theme_gray(base_size = 18)

#Cluster 3 - P.[ST]P
func_list<-filter(as.data.frame(clusters), cluster==3)
func_list<-rownames(func_list)
ann1_temp<-filter(as.data.frame(clusters),cluster==3)
ann1<-data.frame(rownames(ann1_temp),ann1_temp$`.....P.[ST]P......`)
names(ann1)<-c("Description","1-FDR")
func_chart<-read.csv(paste0("Motif_analysis_2/CP_gsb_motif_seqs.txt.....P.[ST]P.......csv"),header=TRUE)
func_chart<-filter(func_chart, Description %in% func_list) 
ann1<-dplyr::full_join(ann1, subset(func_chart, select=c("Description", "GeneRatio","Count")), by="Description")
ann1$GeneRatio<-DOSE::parse_ratio(ann1$GeneRatio)
ann1<-filter(ann1, Count>3)
ann1<-filter(ann1, `1-FDR`>=0.6)
ggplot(ann1, aes(x=GeneRatio, y=Description, colour=`1-FDR`, size=Count))+geom_point()+cs+xlab("Enrichment")+ylab("")+ scale_size_continuous(limits=c(0,60),breaks=c(10,20,30,40,50))+theme_gray(base_size = 18)


#Cluster 5 - R..[ST]
func_list<-filter(as.data.frame(clusters), cluster==5)
func_list<-rownames(func_list)
ann1_temp<-filter(as.data.frame(clusters),cluster==5)
ann1<-data.frame(rownames(ann1_temp),ann1_temp$`....R..[ST].......`)
names(ann1)<-c("Description","1-FDR")
func_chart<-read.csv(paste0("Motif_analysis_2/CP_gsb_motif_seqs.txt....R..[ST]........csv"),header=TRUE)
func_chart<-filter(func_chart, Description %in% func_list) 
ann1<-dplyr::full_join(ann1, subset(func_chart, select=c("Description", "GeneRatio","Count")), by="Description")
ann1$GeneRatio<-DOSE::parse_ratio(ann1$GeneRatio)
ann1<-filter(ann1, Count>3)
ann1<-filter(ann1, `1-FDR`>=0.6)
ggplot(ann1, aes(x=GeneRatio, y=Description, colour=`1-FDR`, size=Count))+geom_point()+cs+xlab("Enrichment")+ylab("")+ scale_size_continuous(limits=c(0,60),breaks=c(10,20,30,40,50))+theme_gray(base_size = 18)



#Cluster 2 - R..[ST]P
func_list<-filter(as.data.frame(clusters), cluster==2)
func_list<-rownames(func_list)
ann1_temp<-filter(as.data.frame(clusters),cluster==2)
ann1<-data.frame(rownames(ann1_temp),ann1_temp$`....R..[ST]P......`)
names(ann1)<-c("Description","1-FDR")
func_chart<-read.csv(paste0("Motif_analysis_2/CP_gsb_motif_seqs.txt....R..[ST]P.......csv"),header=TRUE)
func_chart<-filter(func_chart, Description %in% func_list) 
ann1<-dplyr::full_join(ann1, subset(func_chart, select=c("Description", "GeneRatio","Count")), by="Description")
ann1$GeneRatio<-DOSE::parse_ratio(ann1$GeneRatio)
ann1<-filter(ann1, Count>3)
ann1<-filter(ann1, `1-FDR`>=0.6)
ggplot(ann1, aes(x=GeneRatio, y=Description, colour=`1-FDR`, size=Count))+geom_point()+cs+xlab("Enrichment")+ylab("")+ scale_size_continuous(limits=c(0,60),breaks=c(10,20,30,40,50))+theme_gray(base_size = 18)

