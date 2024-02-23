library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(ggvenn)
library(dplyr)
library(stringr)

#a)bar chart - count non-redundant and redundant phosphopep and sites
FLR_counts<-read.csv("FLRcounts_noA.csv",check.names=FALSE)
FLR_counts<-select(FLR_counts, -c("No. RAW files","Spectral Count"))
FLR_counts<-melt(FLR_counts, id.vars="Dataset", variable.name="Count")
FLR_counts$Group<-ifelse(grepl("Peptidoform",FLR_counts$Count),"Peptidoform","PSM")
FLR_counts$Count<-str_replace(FLR_counts$Count, "Binomial ", "")
FLR_counts$Count<-str_replace(FLR_counts$Count, "pA", "")

b<-ggplot2::ggplot(FLR_counts, aes(fill=Count, y=value, x=Dataset)) + geom_bar(position='dodge', stat='identity')+
  theme(axis.text.x = element_text(angle = 45 , hjust=1))+ggtitle("a)")+facet_grid(factor(Group, levels=c("PSM","Peptidoform")) ~ ., scales="free_y")+
  scale_fill_manual(values=safe_colorblind_palette)+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
b

#b)Gold silver bronze counts per site, with alanines
counts_A<-read.csv("Unique_protein_SNP_reps_binomial_with_no_choice_all.csv")
cat_counts<-subset(counts_A,select=c("PTM_FLR_category","PTM.residue"))
cat_counts<-melt(table(cat_counts))
level_order<-c("Bronze","Silver","Gold")

c<-ggplot(cat_counts,aes(fill=PTM.residue, y=value, x=factor(PTM_FLR_category,level=level_order)))+ 
  geom_bar(position='dodge', stat='identity')+
  geom_text(aes(label = value),size = 3, vjust = -0.5, position = position_dodge(.9))+
  ylab("Count of Sites")+
  xlab("Category")+
  labs(fill="Phosphosite residue")+
  ggtitle("b)")+theme(text = element_text(size=10))+
  scale_fill_manual(values=safe_colorblind_palette[7:11])+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
c

#c) venn diagram
counts<-read.csv("Unique_protein_SNP_reps_binomial_with_no_choice_all.csv")
counts_noA<-filter(counts,PTM.residue!="A")
dataset_counts<-subset(counts_noA,select=c("MSU","RAP_DB","UP"))

dataset_counts$MSU<-ifelse(dataset_counts$MSU!="",TRUE,FALSE)
dataset_counts$RAP_DB<-ifelse(dataset_counts$RAP_DB!="",TRUE,FALSE)
dataset_counts$UP<-ifelse(dataset_counts$UP!="",TRUE,FALSE)

d<-ggplot(dataset_counts) +
  geom_venn(aes(A = `MSU`, B = `RAP_DB`, C = `UP`), 
            show_percentage = FALSE, text_size=3 ,set_name_size=3,stroke_size=0.2, fill_color=colourblind_palette[-1]) +
  coord_fixed() +
  theme_void() +
  ggtitle("c)")
d

prot_counts<-subset(counts_noA,select=c(Protein_count))
prot_counts<-melt(table(prot_counts))
prot_counts$Category<-as.character(prot_counts$Protein_count)
prot_counts$Category[prot_counts$Protein_count>=10] = ">=10"

prot_counts<-aggregate(prot_counts$value, by=list(Category=prot_counts$Category), FUN=sum)

prot_counts$Category<-factor(prot_counts$Category, levels=c("1","2","3","4","5","6","7","8","9",">=10"))

e<-ggplot(data=prot_counts,aes(y=x, x=Category))+ 
  geom_bar(position='dodge', stat='identity')+
  geom_text(aes(label = x),size = 3, vjust = -0.5, position=position_dodge(width=1))+
  ylab("Count of Sites")+
  xlab("Count of Proteins")+
  ggtitle("d)")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

lm <- rbind(c(1,2),
            c(3,4))
p1<-grid.arrange(grobs=list(b,c,d,e), layout_matrix=lm)

