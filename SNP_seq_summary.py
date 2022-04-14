import pandas as pd
import numpy as np
import os
import gzip
from Bio import SeqIO
from bisect import bisect_right

def take_closest(myNumber, myList):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the largest number.
    """
    pos = bisect_right(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber <= myNumber - before:
        return after
    else:
        return before

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
            and len(lst_cols) > 0
            and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
        col:np.repeat(df[col].values, lens)
        for col in idx_cols},
        index=idx)
           .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                      for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
               .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res

flr_filter=0.05
file="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Rice_phosphosite_matrix_updated"+str(flr_filter)+"_w_protein-pos_SNP.csv"
df=pd.read_csv(file, dtype={'Category1count':str,'Category2count':str,'Category3count':str,'Category4count':str,'Stat_file_exist':str,'Fasta_file_exist':str})

df=df.fillna("-")
df['Protein']=df['Protein'].str.split(";")
df['Stat_file_exist']=df['Stat_file_exist'].str.split(";")
df['Fasta_file_exist']=df['Fasta_file_exist'].str.split(";")
df['Category1count']=df['Category1count'].str.split(";")
df['Category2count']=df['Category2count'].str.split(";")
df['Category3count']=df['Category3count'].str.split(";")
df['Category4count']=df['Category4count'].str.split(";")


df = explode(df, ['Protein', 'Stat_file_exist', 'Fasta_file_exist', 'Category1count', 'Category2count', 'Category3count', 'Category4count'], fill_value='')
df[['Protein', 'Start', 'PTM_pos']] = df['Protein'].str.rsplit('_', 2, expand=True)
#df.to_csv(file.replace(".csv","_expand.csv"))
df=df.loc[df['Stat_file_exist']=="1"]
df=df.loc[df['Fasta_file_exist']=="1"]
#protein_list=df['Protein'].unique()

SNP_file_loc="C:/Users/krams/Dropbox/PTMExchange/Rice/SNP/"
SNP_df=pd.DataFrame(columns=['pos', 'aa', 'vargroup', 'count', 'Protein', 'ref_match', 'PTM_pos_temp'])

#takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
closest_list=[]
counter=1
#for p in protein_list:
df=df.reset_index()
for i in range(len(df)):
#for i in range(100):
    print(str(i)+"/"+str(len(df)))
    p=df.loc[i,"Protein"]
    #print(p)
    #print(str(counter)+"/"+str(len(protein_list)))
    counter+=1
    if os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p+"_different_proteins_in_varieties.fasta.gz")==True:
        fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p+"_different_proteins_in_varieties.fasta.gz"
        sub="out_saaps_all_canonical_MSU/"
    elif os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p+"_different_proteins_in_varieties.fasta.gz")==True:
        fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p+"_different_proteins_in_varieties.fasta.gz"
        sub="out_saaps_msu_non_canonical/"
    else:
        fasta=SNP_file_loc+"out_saaps_rapdb/"+p+"_different_proteins_in_varieties.fasta.gz"
        sub="out_saaps_rapdb/"
    with gzip.open(fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if "reference" in record.id:
                seq=record.seq

    PTM_list=df.loc[df["Protein"]==p]['PTM_pos']

    if "LOC" in p:
        SNP_temp_file=SNP_file_loc+sub+p.split(".")[0]+"_stats.csv"
    else:
        SNP_temp_file=SNP_file_loc+sub+p+"_proteins_in_varieties_stats.csv"

    SNP_temp=pd.read_csv(SNP_temp_file)

    SNP_temp['Protein']=p
    SNP_temp['pos']=SNP_temp['pos'].astype(int)

    for s in range(len(SNP_temp)):
        if SNP_temp.loc[s,'aa']==seq[SNP_temp.loc[s,'pos']-1]:
            SNP_temp.loc[s,'ref_match']="Y"
        else:
            SNP_temp.loc[s,'ref_match']="N"
    for PTM in PTM_list:
        SNP_temp[PTM+"_SiteCategory"]=0
        if len(SNP_temp)>0:
            SNP_temp.loc[SNP_temp['pos'].between(int(PTM)-5,int(PTM)+5),PTM+"_SiteCategory"]=4
            SNP_temp.loc[SNP_temp['pos']==int(PTM)+1,PTM+"_SiteCategory"]=2
            SNP_temp.loc[SNP_temp['pos']==int(PTM)-1,PTM+"_SiteCategory"]=3
            SNP_temp.loc[SNP_temp['pos']==int(PTM),PTM+"_SiteCategory"]=1
        #SNP_df=pd.concat([SNP_df,SNP_temp],ignore_index=True,sort=False)
        SNP_temp.to_csv(SNP_temp_file.replace(".csv","_PTM.csv"))

    PTM_query=df.loc[i,'PTM_pos']

    if len(SNP_temp)>0:
        #closest_PTM=takeClosest(int(PTM_query),SNP_temp['pos'].astype(int).unique())
        closest_PTM=take_closest(int(PTM_query),np.sort(SNP_temp['pos'].unique()))
    else:
        closest_PTM=0
    df.loc[i,"PTM residue"]=df.loc[i,"Peptide"].split("-")[0][int(df.loc[i,"Peptide"].split("-")[1])-1]
    #closest_list.append(closest_PTM)
    df.loc[i,'Closest SNP']=closest_PTM
    SNP_res=""
    maf=""
    if closest_PTM!=0:
        SNP_res_temp=SNP_temp.loc[SNP_temp['pos']==closest_PTM]["aa"].tolist()[0]
        SNP_res+=SNP_res_temp+"->"
        SNP_res+=";".join(SNP_temp.loc[(SNP_temp['pos']==closest_PTM) & (SNP_temp['aa']!=SNP_res_temp)]["aa"].unique())
        print(SNP_res)
        f_temp=SNP_res.split(">")[1]
        for f in f_temp.split(";"):
            count=0
            all=0
            with gzip.open(fasta, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.seq[closest_PTM-1]==f:
                        count+=1
                    all+=1
            freq=count/all
            maf+=str(freq)+";"
    df.loc[i,'SNP_res']=SNP_res
    df.loc[i,"Minor Allele Frequency"]=maf[:-1]
#df['Closest SNP']=closest_list
df['PTM_site_category']=0
df=df.loc[df['PTM residue']!="A"]
df.loc[df['Closest SNP'].between(df['PTM_pos'].astype(int)-5,df['PTM_pos'].astype(int)+5),'PTM_site_category']=4
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int)+1,'PTM_site_category']=2
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int)-1,'PTM_site_category']=3
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int),'PTM_site_category']=1
df.to_csv(file.replace(".csv","_expand_per_PTM.csv"))
SNP_df.to_csv(SNP_file_loc+"all_stat_PTMS.csv")

df=pd.read_csv(file.replace(".csv","_expand_per_PTM.csv"))
#Add in oryzabase
oryzabase=pd.read_csv("C:/Users/krams/Dropbox/PTMExchange/Rice/SNP/OryzabaseGeneListAll_20220326010000.csv")
#dict: MSU_ID: CGSNL Gene Symbol; CGSNL Gene Name

oryzabase['Gene']=oryzabase['CGSNL Gene Symbol']+";"+oryzabase["CGSNL Gene Name"]
oryzabase=oryzabase[oryzabase['MSU ID'].notna()]
oryzabase['MSU']=oryzabase['MSU ID'].str.split(".").str[0]
oryzabase_dict=oryzabase.set_index("MSU").to_dict()['Gene']

df['Oryzabase']=df['Protein'].str.split(".").str[0].map(oryzabase_dict)

df.to_csv(file.replace(".csv","_expand_per_PTM_oryzabase_freq.csv"))


