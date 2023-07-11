import pandas as pd
import numpy as np
import os
import gzip
from Bio import SeqIO
from bisect import bisect_right
import re

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
#file="C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Rice_phosphosite_matrix_binomial_rm_no_choice_"+str(flr_filter)+"_w_protein-pos_scores_FLR_SNP.csv"
file="D:/Dropbox\PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/Rice_phosphosite_matrix_binomial_peptidoform_"+str(flr_filter)+"_w_protein-pos_scores_FLR_SNP_reps.csv"
df=pd.read_csv(file, dtype={'Category1count':str,'Category2count':str,'Category3count':str,'Category4count':str,'Stat_file_exist':str,'Fasta_file_exist':str,'Stat_loc':str,'Fasta_loc':str})

df=df.fillna("-")
df['Protein']=df['Protein'].str.split(";")
df['Stat_file_exist']=df['Stat_file_exist'].str.split(";")
df['Fasta_file_exist']=df['Fasta_file_exist'].str.split(";")
df['Stat_loc']=df['Stat_loc'].str.split(";")
df['Fasta_loc']=df['Fasta_loc'].str.split(";")
df['Category1count']=df['Category1count'].str.split(";")
df['Category2count']=df['Category2count'].str.split(";")
df['Category3count']=df['Category3count'].str.split(";")
df['Category4count']=df['Category4count'].str.split(";")
df['Representative_locus']=df['Representative_locus'].str.split(";")


df = explode(df, ['Protein', 'Representative_locus', 'Stat_file_exist', 'Stat_loc', 'Fasta_file_exist', 'Fasta_loc', 'Category1count', 'Category2count', 'Category3count', 'Category4count'], fill_value='')
df[['Protein', 'Start', 'PTM_pos']] = df['Protein'].str.rsplit('_', 2, expand=True)
df.to_csv(file.replace(".csv","_expand.csv"), index=False)
df=df.loc[df['Stat_file_exist']=="1"]
df=df.loc[df['Fasta_file_exist']=="1"]
#protein_list=df['Protein'].unique()

SNP_file_loc="D:/Dropbox\PTMExchange/Rice/SNP/"
SNP_df=pd.DataFrame(columns=['pos', 'aa', 'vargroup', 'count', 'Protein', 'ref_match', 'PTM_pos_temp'])

#takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
closest_list=[]
counter=1
df['PTM residue']=""
df['Closest SNP']=0
df['SNP_res']=""
#for p in protein_list:
df=df.reset_index()
var_list=["trop","temp","admix","japx","subtrop","aus","aro","ind2","indx","ind1B","ind3","ind1A"]
fam_list=["jap","ind","admx","aus"]
df["trop_ref_freq"]=""
df["temp_ref_freq"]=""
df["admix_ref_freq"]=""
df["japx_ref_freq"]=""
df["subtrop_ref_freq"]=""
df["aus_ref_freq"]=""
df["aro_ref_freq"]=""
df["ind2_ref_freq"]=""
df["indx_ref_freq"]=""
df["ind1B_ref_freq"]=""
df["ind3_ref_freq"]=""
df["ind1A_ref_freq"]=""
prot_dict=pd.read_csv(SNP_file_loc+"/3k_pop_details_UTF-8_v3.tsv",sep="\t")

for i in range(len(df)):
#for i in range(10):
    print(str(i)+"/"+str(len(df)))
    if df.loc[i,'Representative_locus']!="NA" and df.loc[i,'Representative_locus']!="-":
        p=df.loc[i,"Representative_locus"]
    else:
        p=df.loc[i,"Protein"]
    counter+=1
    print(p)
    if os.path.exists(SNP_file_loc + "output_SAAVs/output_SAAVs~/output_SAAVs/" + p + "_different_proteins_in_varieties.fasta.gz") == True:
        fasta = SNP_file_loc + "output_SAAVs/output_SAAVs~/output_SAAVs/" + p + "_different_proteins_in_varieties.fasta.gz"
        sub = "output_SAAVs/output_SAAVs~/output_SAAVs/"
    elif os.path.exists(SNP_file_loc+"output_SAAVs_June2023/"+p+"_different_proteins_in_varieties.fasta.gz")==True:
        fasta=SNP_file_loc+"output_SAAVs_June2023/"+p+"_different_proteins_in_varieties.fasta.gz"
        sub="output_SAAVs_June2023/"
    elif os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p+"_proteins_in_varieties.fasta.gz")==True:
        fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p+"_proteins_in_varieties.fasta.gz"
        sub="out_saaps_all_canonical_MSU/"
    elif os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p+"_proteins_in_varieties.fasta.gz")==True:
        fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p+"_proteins_in_varieties.fasta.gz"
        sub="out_saaps_msu_non_canonical/"
    else:
        fasta=SNP_file_loc+"out_saaps_rapdb/"+p+"_proteins_in_varieties.fasta.gz"
        sub="out_saaps_rapdb/"
    with gzip.open(fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if "reference" in record.id:
                seq=record.seq

    PTM_list=df.loc[df["Protein"]==p]['PTM_pos']

    if "LOC" in p:
        if os.path.exists(SNP_file_loc + sub + p + "_stats.csv"):
            SNP_temp_file=SNP_file_loc+sub+p+"_stats.csv"
        else:
            SNP_temp_file=SNP_file_loc+sub+p.split(".")[0]+"_stats.csv"
        df.loc[i,"Protein_pos"]=p.split(".")[0]
    else:
        SNP_temp_file=SNP_file_loc+sub+p+"_proteins_in_varieties_stats.csv"
        df.loc[i,"Protein_pos"]=p.split("-")[0]

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
        #print(SNP_temp_file)
        SNP_temp.to_csv(SNP_temp_file.replace(".csv","_PTM.csv"), index=False)

    PTM_query=df.loc[i,'PTM_pos']

    if len(SNP_temp)>0:
        #closest_PTM=takeClosest(int(PTM_query),SNP_temp['pos'].astype(int).unique())
        closest_PTM=take_closest(int(PTM_query),np.sort(SNP_temp['pos'].unique()))
    else:
        closest_PTM=0
    #print(df.loc[i,"Peptide"])
    #peptide_no_mod=re.sub("([\(\[]).*?([\)\]])","",df.loc[i,"Peptide"].rsplit("-",1)[0])
    peptide_no_mod=df.loc[i,"Peptide"]
    #(peptide_no_mod)
    df.loc[i,"PTM residue"]=peptide_no_mod[int(df.loc[i,"Peptide_mod"].rsplit("-",1)[1])-1]
    #print(df.loc[i,"PTM residue"])
    #closest_list.append(closest_PTM)
    df.loc[i,'Closest SNP']=closest_PTM
    SNP_res=""
    maf=""
    if closest_PTM!=0:
        SNP_res_temp=SNP_temp.loc[SNP_temp['pos']==closest_PTM]["aa"].tolist()[0]
        SNP_res+=SNP_res_temp+"->"
        SNP_res+=";".join(SNP_temp.loc[(SNP_temp['pos']==closest_PTM) & (SNP_temp['aa']!=SNP_res_temp)]["aa"].unique())
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
    if SNP_res!="":
        for v in var_list:
            prot_dict_temp=prot_dict.loc[prot_dict["SUBPOPULATION"]==v]
            prot_list=prot_dict_temp["ASSAY ID"].tolist()
            count=0
            all=0
            f=SNP_res.split("-")[0]
            with gzip.open(fasta, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if "LOC" not in record.id:
                        query=record.id.split("_")[1]
                    else:
                        query=record.id.split("_",2)[2]
                    if query in prot_list:
                        if record.seq[closest_PTM-1]==f:
                            count+=1
                        all+=1
            if all!=0:
                df.loc[i,v+"_ref_freq"]=count/all
        for f in fam_list:
            df.loc[i,f+"_family_count"]=0
            if f=="jap":
                fam=["trop","temp","japx","subtrop","aro"]
            elif f=="ind":
                fam=["ind2","indx","ind1B","ind3","ind1A"]
            elif f=="admx":
                fam=["admix"]
            else:
                fam=["aus"]
            tot=0
            for f2 in fam:
                if df.loc[i,f2+"_ref_freq"]!="":
                    tot+= float(df.loc[i,f2+"_ref_freq"])
            df.loc[i,f+"_family_ref_freq"]=tot/len(fam)
            df.loc[i,f+"_family_count"]=len(fam)

#diff ["jap","ind","admx","aus"]
df["Jap_Ind_Diff"]=abs(df["jap_family_ref_freq"]-df["ind_family_ref_freq"])
df["Jap_Admx_diff"]=abs(df["jap_family_ref_freq"]-df["admx_family_ref_freq"])
df["Jap_Aus_diff"]=abs(df["jap_family_ref_freq"]-df["aus_family_ref_freq"])
df["Ind_Admx_diff"]=abs(df["ind_family_ref_freq"]-df["admx_family_ref_freq"])
df["Ind_Aus_diff"]=abs(df["ind_family_ref_freq"]-df["aus_family_ref_freq"])
df["Admx_Aus_diff"]=abs(df["admx_family_ref_freq"]-df["aus_family_ref_freq"])
df["Total_diff"]=df["Jap_Ind_Diff"]+df["Jap_Admx_diff"]+df["Jap_Aus_diff"]+df["Ind_Admx_diff"]+df["Ind_Aus_diff"]+df["Admx_Aus_diff"]
#df['Closest SNP']=closest_list
df['PTM_site_category']=0
#df=df.loc[df['PTM residue']!="A"]
df.loc[df['Closest SNP'].between(df['PTM_pos'].astype(int)-5,df['PTM_pos'].astype(int)+5),'PTM_site_category']=4
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int)+1,'PTM_site_category']=2
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int)-1,'PTM_site_category']=3
df.loc[df['Closest SNP']==df['PTM_pos'].astype(int),'PTM_site_category']=1
df.to_csv(file.replace(".csv","_expand_per_PTM.csv"), index=False)
#SNP_df.to_csv(SNP_file_loc+"all_stat_PTMS.csv", index=False)

df['Protein_pos']=df['Protein_pos']+"-"+df['PTM_pos']
df['Protein_count'] = df.groupby('Protein_pos')['Protein_pos'].transform('size')
df['Unique']=df['Protein_count']
df.loc[df['Unique'] > 1, 'Unique'] = 0


#Add in oryzabase
oryzabase=pd.read_csv("D:/Dropbox\PTMExchange/Rice/SNP/OryzabaseGeneListAll_20220326010000.csv")
#dict: MSU_ID: CGSNL Gene Symbol; CGSNL Gene Name
oryzabase['Gene']=oryzabase['CGSNL Gene Symbol']+";"+oryzabase["CGSNL Gene Name"]
oryzabase=oryzabase[oryzabase['MSU ID'].notna()]
oryzabase['MSU']=oryzabase['MSU ID'].str.split(".").str[0]
oryzabase_dict=oryzabase.set_index("MSU").to_dict()['Gene']
#df['Oryzabase']=df['Protein'].str.split(".").str[0].map(oryzabase_dict)
df['Oryzabase']=df['MSU'].str.split(".").str[0].map(oryzabase_dict)

#Add in RAP-DB genes
RAP_genes=pd.read_csv("D:/Dropbox\PTMExchange/Rice/SNP/RAP_DB_IRGSP-1.0_representative_annotation_2022-03-11/IRGSP-1.0_representative_annotation_2022-03-11.tsv",sep="\t")
#dict: RAP_ID: CGSNL Gene Symbol; CGSNL Gene Name; Desription
RAP_genes['RAP-DB Gene Symbol Synonym(s)'] = RAP_genes['RAP-DB Gene Symbol Synonym(s)'].fillna("-")
RAP_genes['RAP-DB Gene Name Synonym(s)'] = RAP_genes['RAP-DB Gene Name Synonym(s)'].fillna("-")
RAP_genes['Gene']=RAP_genes['RAP-DB Gene Symbol Synonym(s)']+";"+RAP_genes["RAP-DB Gene Name Synonym(s)"]+";"+RAP_genes["Description"]
#RAP_genes=RAP_genes[RAP_genes['Transcript_ID'].notna()]
RAP_genes=RAP_genes.set_index(RAP_genes["Transcript_ID"])
RAP_genes_dict=RAP_genes.to_dict()['Gene']
df['RAP_genes']=df['RAP_DB'].str.split("_").str[0].map(RAP_genes_dict)

#Add in MSU annotation
MSU_genes=pd.read_csv("D:/Dropbox\PTMExchange/Rice/SNP/MSU_pfam.tsv",sep="\t")
MSU_genes['hmm_acc'] = MSU_genes['hmm_acc'].fillna("-")
MSU_genes['hmm_name'] = MSU_genes['hmm_name'].fillna("-")
MSU_genes['Gene']=MSU_genes['hmm_acc']+";"+MSU_genes["hmm_name"]+";"+MSU_genes["hmm_type"]
MSU_genes=MSU_genes.set_index(MSU_genes["model"])
MSU_genes_dict=MSU_genes.to_dict()['Gene']
df['MSU_genes']="LOC_"+df['MSU'].str.split("_").str[1]
df['MSU_genes']=df['MSU_genes'].map(MSU_genes_dict)

#cols=['Oryzabase','RAP_genes','MSU_genes']
#df[cols]=df[cols].fillna(0)

df['Annotated?']=False
#df.loc[df['MSU_genes'].str.contains('[a-zA-Z]',na=False),'Annotated?']=True
#df.loc[df['Oryzabase'].str.contains('[a-zA-Z]',na=False),'Annotated?']=True
#df.loc[df['RAP_genes'].str.contains('[a-zA-Z]',na=False),'Annotated?']=True
df.loc[(df['MSU_genes'].str.contains('[a-zA-Z]',na=False))|(df['RAP_genes'].str.contains('[a-zA-Z]',na=False))|(df['Oryzabase'].str.contains('[a-zA-Z]',na=False)),'Annotated?']=True

df.to_csv(file.replace(".csv","_expand_per_PTM_oryzabase_RAP-DB_MSU_freq_var.csv"), index=False)

score_cols = [col for col in df.columns if '_score' in col]
df[score_cols]=df[score_cols].replace("-",0)
df[score_cols] = df[score_cols].apply(pd.to_numeric)
df['Average_score']=df[score_cols].sum(axis=1)/len(score_cols)
df=df.sort_values(by='Average_score', ascending=False)
df=df.drop_duplicates(subset=['Protein_pos'], keep='first')

#df.to_csv("C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Unique_protein_SNP_binomial.csv", index=False)
df.to_csv("D:/Dropbox\PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/Unique_protein_SNP_reps_binomial_with_no_choice.csv", index=False)