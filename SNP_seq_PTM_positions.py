import pandas as pd
import os
import gzip
from Bio import SeqIO


flr_filter=0.05
file="C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/Rice_phosphosite_matrix_updated"+str(flr_filter)+"_w_protein-pos.csv"
#file="C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/Rice_phosphosite_matrix_updated"+str(flr_filter)+".csv"
df=pd.read_csv(file)

SNP_file_loc="C:/Users/krams/Dropbox/PTMExchange/Rice/SNP/"

stat_file_list=[]
fasta_list=[]

SNP_df=pd.DataFrame(columns=['pos', 'aa', 'vargroup', 'count', 'Protein', 'ref_match', 'PTM_pos_temp', 'PTM_site_category'])
#category1: ST switch, lose PTM site (eg. S -> A)
#category2: motif site switch +1 position, disrupt phospho motif (eg. S*)
#category3: motif site switch -1 position, disrupt phospho motif (eg. *S)
#category4: polymorphism with +/-5aa pf site, potentially motifies site

for i in range(len(df)):
    #11:55
    print(i, len(df))
#for i in range(100):
    if pd.notna(df.loc[i,"MSU"]):
        protein=df.loc[i,"MSU"]
    else:
        protein=df.loc[i,"RAP_DB"]
    if pd.notna(protein):
        for p in protein.split(";"):
            PTMpos=int(p.rsplit("_",1)[1])
            if "LOC" in p:
                fasta="NOT_FOUND"
                if os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
                elif os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
                else:
                    if p.rsplit('_', 1)[0] not in fasta_list:
                        fasta_list.append(p.rsplit('_', 2)[0])
                if fasta!="NOT_FOUND":
                    with gzip.open(fasta, "rt") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            if "reference" in record.id:
                                seq=record.seq
                else:
                    seq="NA"

                try:
                    SNP_temp=pd.read_csv(SNP_file_loc+"Stats_files_all/"+p.split(".")[0]+"_stats.csv")
                    SNP_temp['Protein']=p.split(".")[0]
                    SNP_temp['pos']=SNP_temp['pos'].astype(int)
                    if seq=="NA":
                        SNP_temp['ref_match']="NA"
                    else:
                        SNP_temp['ref_match']="N"
                        #SNP_temp.loc[SNP_temp['aa']==seq[SNP_temp['pos']-1],'ref_match']="Y"
                    SNP_temp['PTM_pos_temp']=PTMpos
                    SNP_temp['PTM_site_category']=0
                    SNP_temp.loc[SNP_temp['pos'].between(PTMpos-5,PTMpos+5),'PTM_site_category']=4
                    SNP_temp.loc[SNP_temp['pos']==PTMpos+1,'PTM_site_category']=2
                    SNP_temp.loc[SNP_temp['pos']==PTMpos-1,'PTM_site_category']=3
                    SNP_temp.loc[SNP_temp['pos']==PTMpos,'PTM_site_category']=1

                    for s in range(len(SNP_temp)):
                        if SNP_temp.loc[s,'aa']==seq[SNP_temp.loc[s,'pos']-1]:
                            SNP_temp.loc[s,'ref_match']="Y"
                        else:
                            SNP_temp.loc[s,'ref_match']="N"

                    SNP_df=pd.concat([SNP_df,SNP_temp],ignore_index=True,sort=False)
                except:
                    if p.rsplit('_', 1)[0] not in stat_file_list:
                        stat_file_list.append(p.rsplit('_', 1)[0])
            elif "Os" in p:
                fasta="NOT_FOUND"
                if os.path.exists(SNP_file_loc+"out_saaps_rapdb/"+p.split("_")[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_rapdb/"+p.split("_")[0]+"_different_proteins_in_varieties.fasta.gz"
                else:
                    if p.rsplit('_', 1)[0] not in fasta_list:
                        fasta_list.append(p.split('_')[0])
                if fasta!="NOT_FOUND":
                    with gzip.open(fasta, "rt") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            if "reference" in record.id:
                                seq=record.seq
                else:
                    seq="NA"

                try:
                    SNP_temp=pd.read_csv(SNP_file_loc+"Stats_files_all/"+p.split("_")[0]+"_proteins_in_varieties_stats.csv")
                    SNP_temp['Protein']=p.split(".")[0]
                    SNP_temp['pos']=SNP_temp['pos'].astype(int)
                    if seq=="NA":
                        SNP_temp['ref_match']="NA"
                    else:
                        SNP_temp['ref_match']="N"
                        #SNP_temp.loc[SNP_temp['aa']==seq[SNP_temp['pos']-1],'ref_match']="Y"
                    SNP_temp['PTM_pos_temp']=PTMpos
                    SNP_temp['PTM_site_category']=0
                    SNP_temp.loc[SNP_temp['pos'].between(PTMpos-5,PTMpos+5),'PTM_site_category']=4
                    SNP_temp.loc[SNP_temp['pos']==PTMpos+1,'PTM_site_category']=2
                    SNP_temp.loc[SNP_temp['pos']==PTMpos-1,'PTM_site_category']=3
                    SNP_temp.loc[SNP_temp['pos']==PTMpos,'PTM_site_category']=1

                    for s in range(len(SNP_temp)):
                        if SNP_temp.loc[s,'aa']==seq[SNP_temp.loc[s,'pos']-1]:
                            SNP_temp.loc[s,'ref_match']="Y"
                        else:
                            SNP_temp.loc[s,'ref_match']="N"
                            print(seq,SNP_temp.loc[s,'pos'],SNP_temp.loc[s,'aa'])

                    SNP_df=pd.concat([SNP_df,SNP_temp],ignore_index=True,sort=False)
                except:
                    if p.rsplit('_', 1)[0] not in stat_file_list:
                        stat_file_list.append(p.rsplit('_', 2)[0])
            else:
                print("Check protein: "+p)
print(stat_file_list)
print(fasta_list)
SNP_df.to_csv(SNP_file_loc+"all_stat_PTMS.csv")


