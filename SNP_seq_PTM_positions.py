import pandas as pd
import os
import gzip
from Bio import SeqIO


flr_filter=0.05
#file="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Rice_phosphosite_matrix_binomial_rm_no_choice_"+str(flr_filter)+"_w_protein-pos_scores_FLR.csv"
file="D:/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/Rice_phosphosite_matrix_binomial_peptidoform_"+str(flr_filter)+"_w_protein-pos_scores_FLR.csv"
#file="C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/Rice_phosphosite_matrix_updated"+str(flr_filter)+".csv"
df=pd.read_csv(file)

SNP_file_loc="D:/Dropbox/PTMExchange/Rice/SNP/"

transcripts = pd.read_csv("D:/Dropbox/PTMExchange/Rice/SNP/all.locus_brief_info.7.0", sep="\t") #tsv file
transcripts = transcripts.loc[transcripts['is_representative']=="Y"]
transcripts_dict = dict(zip(transcripts.locus, transcripts.model))

stat_file_list=[]
fasta_list_all=[]

SNP_df=pd.DataFrame(columns=['pos', 'aa', 'vargroup', 'count', 'Protein', 'ref_match', 'PTM_pos_temp', 'PTM_site_category'])
#category1: ST switch, lose PTM site (eg. S -> A)
#category2: motif site switch +1 position, disrupt phospho motif (eg. S*)
#category3: motif site switch -1 position, disrupt phospho motif (eg. *S)
#category4: polymorphism with +/-5aa pf site, potentially motifies site

df['Protein']=0
df['Stat_file_exist']=0
df['Stat_loc']=0
df['Fasta_file_exist']=0
df['Fasta_loc']=0
df['Category1count']=0
df['Category2count']=0
df['Category3count']=0
df['Category4count']=0
for i in range(len(df)):
    protein=""
    print(i, len(df))
#for i in range(100):
    cat1_list=""
    cat2_list=""
    cat3_list=""
    cat4_list=""
    stat_list=""
    fasta_list=""
    stat_loc_list=""
    fasta_loc_list=""

    # if not mapped to MSU - use RAP-DB
    #if pd.notna(df.loc[i,"MSU"]):
    #    protein=df.loc[i,"MSU"]
    #else:
    #    protein=df.loc[i,"RAP_DB"]

    # Only use MSU mapped
    if pd.notna(df.loc[i,"MSU"]):
        protein=df.loc[i,"MSU"]
    df.loc[i,'Protein']=protein

    reps=""
    if pd.notna(protein):
        for p in protein.split(";"):
            if "." in p and "Chr" not in p:
                rep=transcripts_dict[p.split(".")[0]]
                #reps+=rep+";"
            #else:
                #reps+="NA;"
            if p!="":
                PTMpos=int(p.rsplit("_",1)[1])
            if "LOC" in p:
                fasta="NOT_FOUND"
                if os.path.exists(SNP_file_loc+"output_SAAVs/output_SAAVs~/output_SAAVs/"+rep+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"output_SAAVs/output_SAAVs~/output_SAAVs/"+rep+"_different_proteins_in_varieties.fasta.gz"
                    sub="output_SAAVs/output_SAAVs~/output_SAAVs/"+rep
                    reps += rep + ";"
                elif os.path.exists(SNP_file_loc+"output_SAAVs_June2023/"+rep+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"output_SAAVs_June2023/"+rep+"_different_proteins_in_varieties.fasta.gz"
                    sub="output_SAAVs_June2023/"+rep
                    reps += rep + ";"
                elif os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
                    sub="out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]
                    reps += "NA;"
                elif os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+rep+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+rep+"_different_proteins_in_varieties.fasta.gz"
                    sub="out_saaps_all_canonical_MSU/"+rep
                    reps += rep+";"
                elif os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
                    sub="out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]
                    reps += "NA;"
                else:
                    sub=""
                    if p.rsplit('_', 1)[0] not in fasta_list:
                        fasta_list_all.append(p.rsplit('_', 1)[0])
                        reps += rep+";"

                if fasta!="NOT_FOUND":
                    fasta_list+="1;"
                    fasta_loc_list+=fasta+";"
                    with gzip.open(fasta, "rt") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            if "reference" in record.id:
                                seq=record.seq

                else:
                    seq="NA"
                    fasta_list+="0;"
                    fasta_loc_list+="-;"

                try:
                    if os.path.exists(SNP_file_loc+sub+"_stats.csv"):
                        SNP_temp = pd.read_csv(SNP_file_loc + sub + "_stats.csv")
                    else:
                        SNP_temp=pd.read_csv(SNP_file_loc+sub.rsplit(".",1)[0]+"_stats.csv")
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
                        try:
                            if SNP_temp.loc[s,'aa']==seq[SNP_temp.loc[s,'pos']-1]:
                                SNP_temp.loc[s,'ref_match']="Y"
                            else:
                                SNP_temp.loc[s,'ref_match']="N"
                        except:
                            SNP_temp.loc[s,'ref_match']="N"

                    SNP_df=pd.concat([SNP_df,SNP_temp],ignore_index=True,sort=False)
                    stat_list+="1;"
                    stat_loc_list+=SNP_file_loc+sub+"_stats.csv;"
                    SNP_temp=SNP_temp.drop_duplicates(subset = ['pos', 'PTM_site_category'],keep = 'last').reset_index(drop = True)
                    cat1_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==1]))+";"
                    cat2_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==2]))+";"
                    cat3_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==3]))+";"
                    cat4_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==4]))+";"
                except:
                    if p.rsplit('_', 1)[0] not in stat_file_list:
                        stat_file_list.append(p.rsplit('_', 1)[0])
                    stat_list+="0;"
                    stat_loc_list+="-;"
                    cat1_list+="0;"
                    cat2_list+="0;"
                    cat3_list+="0;"
                    cat4_list+="0;"
            elif "Os" in p:
                reps += "NA;"
                fasta="NOT_FOUND"
                if os.path.exists(SNP_file_loc+"out_saaps_rapdb/"+p.split("_")[0]+"_different_proteins_in_varieties.fasta.gz")==True:
                    fasta=SNP_file_loc+"out_saaps_rapdb/"+p.split("_")[0]+"_different_proteins_in_varieties.fasta.gz"
                    sub="out_saaps_rapdb/"
                    fasta_list+="1;"
                    fasta_loc_list+=fasta+";"
                else:
                    if p.rsplit('_', 1)[0] not in fasta_list:
                        fasta_list_all.append(p.split('_')[0])
                if fasta!="NOT_FOUND":
                    with gzip.open(fasta, "rt") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            if "reference" in record.id:
                                seq=record.seq
                else:
                    seq="NA"
                    fasta_list+="0;"
                    fasta_loc_list+="-;"

                try:
                    SNP_temp=pd.read_csv(SNP_file_loc+sub+p.split("_")[0]+"_proteins_in_varieties_stats.csv")
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
                    stat_list+="1;"
                    stat_loc_list+=SNP_file_loc+sub+p.split("_")[0]+"_proteins_in_varieties_stats.csv;"
                    SNP_temp=SNP_temp.drop_duplicates(subset = ['pos', 'PTM_site_category'],keep = 'last').reset_index(drop = True)
                    cat1_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==1]))+";"
                    cat2_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==2]))+";"
                    cat3_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==3]))+";"
                    cat4_list+=str(len(SNP_temp.loc[SNP_temp['PTM_site_category']==4]))+";"
                except:
                    if p.rsplit('_', 1)[0] not in stat_file_list:
                        stat_file_list.append(p.rsplit('_', 2)[0])
                    stat_list+="0;"
                    stat_loc_list+="-;"
                    cat1_list+="0;"
                    cat2_list+="0;"
                    cat3_list+="0;"
                    cat4_list+="0;"
            else:
                print("Check protein: "+p)
                fasta_list+="0;"
                fasta_loc_list+="-;"
                stat_list+="0;"
                stat_loc_list+="-;"
                cat1_list+="0;"
                cat2_list+="0;"
                cat3_list+="0;"
                cat4_list+="0;"
                reps += "NA;"

    df.loc[i,'Representative_locus']=reps[:-1]
    df.loc[i,'Fasta_file_exist']=fasta_list[:-1]
    df.loc[i,'Fasta_loc']=fasta_loc_list[:-1]
    df.loc[i,'Stat_file_exist']=stat_list[:-1]
    df.loc[i,'Stat_loc']=stat_loc_list[:-1]
    df.loc[i,'Category1count']=cat1_list[:-1]
    df.loc[i,'Category2count']=cat2_list[:-1]
    df.loc[i,'Category3count']=cat3_list[:-1]
    df.loc[i,'Category4count']=cat4_list[:-1]

#print(stat_file_list)
#print(fasta_list_all)
df.to_csv(file.replace(".csv","_SNP_reps.csv"))
SNP_df.to_csv(SNP_file_loc+"all_stat_PTMS.csv")


