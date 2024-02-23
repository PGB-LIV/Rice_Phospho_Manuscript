# Analysis of Rice phosphoprotome build 

Files for re-analysis of rice phosphoproteomics datasets located in PRIDE [PXD046118](https://www.ebi.ac.uk/pride/archive/projects/PXD046188) using pAla FLR pipeline

https://github.com/PGB-LIV/mzidFLR

## Genomic site data with SNP annotation

To generate supplementary file *Gemomic_site_data_w_SNP_annotation.csv*:

Uses outputs from following repros:
1.  https://github.com/andrewrobertjones/rice_phospho_ptm_paper
2.  https://github.com/CBFLivUni/mapping_and_enrichment

Requirements:
+ [RAP_DB sequence fasta](https://rapdb.dna.affrc.go.jp/download/irgsp1.html)
+ [MSU sequence fasta](http://rice.uga.edu/downloads_gad.shtml)
+ [Uniprot sequence fasta](https://www.uniprot.org/uniprotkb?query=oryza+sativa&facets=model_organism%3A39947)

### Python scripts:
![Workflow image](https://github.com/kramsbottom/Rice_scripts/assets/57440286/441379a6-6e30-4d34-8254-fcea7590db35)

**Motif_seqs.py** 

Generates sequences used for motif analysis 

Requirements:
+ Output from *SNP_seq_summary.py*
+ Database used for re-analysis results [*ie. Osativa_super_annotation_union_noIC4R_v2_cRAP.fasta*](https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/11/PXD046188/Osativa_super_annotation_union_noIC4R_v2_cRAP_decoy.fasta)

### R scripts:
To generate manuscript figures:
- Figure 1: Overview of phosphosite counts.
- Figure 3: Motif analysis.
- Figure 4: SAAV analysis.

Located in Figs subfolder
  
