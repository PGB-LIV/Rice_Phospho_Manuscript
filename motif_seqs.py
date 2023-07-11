import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os

all_seq_fg=[]
all_seq_prot_fg=[]
all_motif_bg=[]
all_motif_prot_bg=[]
all_phospho_bg=[]

phospho_sites="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/Rice_phosphosite_matrix_binomial_peptidoform_0.05_w_protein-pos_scores_FLR_SNP_reps_expand.csv"
database="C:/Users/krams/Dropbox/PTMExchange/Rice/database/protein_set_comparisons/protein_set_comparisons/Osativa_super_annotation_union_noIC4R_v2_cRAP.fasta"

df = pd.read_csv(phospho_sites)
#list for all proteins, remove duplicated protein phosphosites
df=df.dropna(subset = ['PTM_pos'], inplace = False)
df=df.loc[df['UP']!="-"]
df['PTM_pos']=df['PTM_pos'].astype(int)
PTM_pos_list=df['PTM_pos'].tolist()
protein_list=df['Protein'].to_list()
MSU_list=df["MSU"].to_list()
RAP_DB_list=df["RAP_DB"].to_list()
peptide_list=df['Peptide'].to_list()
category_list=df['PTM_FLR_category'].to_list()

seq_dict=SeqIO.to_dict(SeqIO.parse(database,"fasta",alphabet=IUPAC.extended_protein))
seq_list=[]
seq_protein=[]
bronze_seq_list=[]
bronze_seq_protein=[]
silver_seq_list=[]
silver_seq_protein=[]
gold_seq_list=[]
gold_seq_protein=[]

#All STY phosphosites +/- 7 - 15mer peptides (motif foreground)
for a,b,c,d,e,f in zip(PTM_pos_list,protein_list,peptide_list,category_list,MSU_list,RAP_DB_list):
    e=e.rsplit("_",2)[0]
    f=f.split("_")[0]
    if b in seq_dict:
        record=seq_dict[b]
    elif "sp|"+b in seq_dict:
        record=seq_dict["sp|"+b]
    elif "tr|"+b in seq_dict:
        record=seq_dict["tr|"+b]
    elif e in seq_dict:
        record=seq_dict[e]
    elif f in seq_dict:
        record=seq_dict[f]
    else:
        print(b,e,f)
    seq_temp=record.seq

    if a+9>len(seq_temp):
        while a+9>len(seq_temp):
            seq_temp+="_"
    if a<8:
        while a<8:
            seq_temp="_"+seq_temp
            a+=1
    seq_list.append(seq_temp[a-8:a+7])
    all_seq_fg.append(seq_temp[a-8:a+7])

    seq_protein.append(b)

    if d=="Bronze":
        bronze_seq_list.append(seq_temp[a-8:a+7])
        bronze_seq_protein.append(b)
    if d=="Silver":
        silver_seq_list.append(seq_temp[a-8:a+7])
        silver_seq_protein.append(b)
    if d=="Gold":
        gold_seq_list.append(seq_temp[a-8:a+7])
        gold_seq_protein.append(b)

#save as list of sequences - look up Uniprot (DAVID foreground)
UP_dict=dict(zip(df.Protein, df.UP))
all_seq_list=[seq_list,bronze_seq_list,silver_seq_list,gold_seq_list]
all_proteins=[seq_protein,bronze_seq_protein,silver_seq_protein,gold_seq_protein]
names=["motif_seqs.txt","bronze_motif_seqs.txt","silver_motif_seqs.txt","gold_motif_seqs.txt"]
for i,j,k in zip(all_seq_list,all_proteins,names):
    df2 = pd.DataFrame(list(zip(i,j)), columns=['Sequences','Proteins'])
    updated_protein_list=[]
    for a in range(len(df2)):
        protein=df2.loc[a,'Proteins']
        #if protein isn't uniprot accession ("|")
        protein_2=UP_dict[protein].split("|")[1]
        print(protein,protein_2)
        updated_protein_list.append(protein_2)
        if k=="motif_seqs.txt":
            all_seq_prot_fg.append(protein_2)

    df2['Protein']=updated_protein_list
    df2=df2.drop_duplicates()
    df2.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/"+k,index=False, header=False)


#All STY sites +/- 7 - 15mer (motif background)
df3 = pd.read_csv(phospho_sites)
df3=df3.dropna(subset = ['PTM_pos'], inplace = False)
df3=df3.loc[df3['UP']!="-"]
df3['PTM_pos']=df3['PTM_pos'].astype(int)
df3['Start']=df3['Start'].astype(int)
protein_list2=df3['Protein'].to_list()
peptide_list2=df3['Peptide'].to_list()
protein_start_list=df3['Start'].to_list()
background_list=[]
background_protein=[]
for a,b,c in zip(protein_list2,peptide_list2,protein_start_list):
    if a in seq_dict:
        record=seq_dict[a]
    elif "sp|"+a in seq_dict:
        record=seq_dict["sp|"+a]
    elif "tr|"+a in seq_dict:
        record=seq_dict["tr|"+a]
    else:
        print(a)
    seq_temp=record.seq
    print(seq_temp)
    pos=c
    print(pos)
    for i in b.split("-")[0]:
        if i=="S" or i=="T" or i=="Y":
            if pos+9>len(seq_temp):
                while pos+9>len(seq_temp):
                    seq_temp+="_"
            if pos<8:
                while pos<8:
                    seq_temp="_"+seq_temp
                    pos+=1
            background_list.append(seq_temp[pos-8:pos+7])
            all_motif_bg.append(seq_temp[pos-8:pos+7])
            background_protein.append(a)
            all_motif_prot_bg.append(a)
        pos+=1

UP_dict=dict(zip(df.Protein, df.UP))
all_seq_list=[seq_list,bronze_seq_list,silver_seq_list,gold_seq_list]
all_proteins=[seq_protein,bronze_seq_protein,silver_seq_protein,gold_seq_protein]
names=["motif_seqs.txt","bronze_motif_seqs.txt","silver_motif_seqs.txt","gold_motif_seqs.txt"]
print(background_protein)

df2 = pd.DataFrame(list(zip(background_list,background_protein)), columns=['Sequences','Proteins'])
updated_protein_list=[]
for a in range(len(df2)):
    protein=df2.loc[a,'Proteins']
    print(protein)
    #if protein isn't uniprot accession ("|")
    protein_2=UP_dict[protein].split("|")[1]
    print(protein,protein_2)
    updated_protein_list.append(protein_2)

print(updated_protein_list)

df4 = pd.DataFrame(list(zip(background_list,background_protein)), columns=['Sequences','Proteins'])
df4['Protein2']=updated_protein_list
print(df4)
df4=df4.drop_duplicates()
df4.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/motif_background.txt",index=False, header=False)

#All phosphopeptide proteins - look up Uniprot (DAVID background)
updated_protein_list=[]
for protein in (protein_list):
    #if protein isn't uniprot accession ("|")
    protein_2=UP_dict[protein].split("|")[1]
    updated_protein_list.append(protein_2)
    all_phospho_bg.append(protein_2)

df5=pd.DataFrame(list(updated_protein_list),columns=['Proteins'])
df5=df5.drop_duplicates()
df5.to_csv("All_phospoproteins_background.txt",index=False, header=False)

df6=pd.DataFrame(list(zip(all_seq_fg,all_seq_prot_fg)), columns=['Sequences','Proteins'])
df6=df6.drop_duplicates()
df6.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/motif_seqs.txt",index=False, header=False)
#df7 = pd.DataFrame(list(zip(all_motif_bg,all_motif_prot_bg)), columns=['Sequences','Proteins'])
#df7.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/All_datasets/motif_background.txt",index=False, header=False)
df8=pd.DataFrame(list(all_phospho_bg),columns=['Proteins'])
df8=df8.drop_duplicates()
df8.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated//05.05.23/All_datasets/All_phospoproteins_background.txt",index=False, header=False)