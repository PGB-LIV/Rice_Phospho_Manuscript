import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import os
import gzip

def extract_peptides(input):
    dict = {}
    for seq_record in SeqIO.parse(input, "fasta"):
        seq = str(seq_record.seq).strip("*")
        protein = seq_record.id
        # cut at K or R followed by P
        pep_temp = re.sub(r'(?<=[RK])(?=[^P])', '\n', seq, 0, re.DOTALL)
        pep = (pep_temp.split())
        protein_length = 1
        pos = []
        peptides = []
        positions = []

        # start position of peptides
        for p in pep:
            start = protein_length
            protein_length += len(p)
            pos.append(start)
        # add peptides and start pos to dictionary
        for i in zip(pep, pos):
            peptide = (i[0])
            peptides.append(peptide)
            position = (i[1])
            positions.append(position)
            if len(peptide) >= 5:
                if peptide not in dict.keys():
                    dict[peptide] = protein + "_" + str(position)
                else:
                    temp = str(dict[peptide])
                    if protein + "_" + str(position) not in temp:
                        dict[peptide] = temp + ";" + protein + "_" + str(position)
        # add missed cleavages - concat (max 4 per peptide)
        count = 0
        for i in peptides:
            pos2 = count
            count += 1
            if pos2 <= (len(peptides) - 5):
                missed4 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2] + peptides[pos2 + 3] + peptides[
                    pos2 + 4])
                missed3 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2] + peptides[pos2 + 3])
                missed2 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2])
                missed = (peptides[pos2] + peptides[pos2 + 1])
                pos_missed = (positions[pos2])
                if missed4 not in dict.keys():
                    dict[missed4] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed4])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed4] = temp + ";" + protein + "_" + str(pos_missed)
                if missed3 not in dict.keys():
                    dict[missed3] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed3])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed3] = temp + ";" + protein + "_" + str(pos_missed)
                if missed2 not in dict.keys():
                    dict[missed2] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed2])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed2] = temp + ";" + protein + "_" + str(pos_missed)
                if missed not in dict.keys():
                    dict[missed] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed] = temp + ";" + protein + "_" + str(pos_missed)
            if pos2 <= (len(peptides) - 4):
                missed3 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2] + peptides[pos2 + 3])
                missed2 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2])
                missed = (peptides[pos2] + peptides[pos2 + 1])
                pos_missed = (positions[pos2])
                if missed3 not in dict.keys():
                    dict[missed3] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed3])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed3] = temp + ";" + protein + "_" + str(pos_missed)
                if missed2 not in dict.keys():
                    dict[missed2] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed2])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed2] = temp + ";" + protein + "_" + str(pos_missed)
                if missed not in dict.keys():
                    dict[missed] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed] = temp + ";" + protein + "_" + str(pos_missed)
            if pos2 <= (len(peptides) - 3):
                missed2 = (peptides[pos2] + peptides[pos2 + 1] + peptides[pos2 + 2])
                missed = (peptides[pos2] + peptides[pos2 + 1])
                pos_missed = (positions[pos2])
                if missed2 not in dict.keys():
                    dict[missed2] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed2])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed2] = temp + ";" + protein + "_" + str(pos_missed)
                if missed not in dict.keys():
                    dict[missed] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed] = temp + ";" + protein + "_" + str(pos_missed)
            if pos2 <= (len(peptides) - 2):
                pos_missed = (positions[pos2])
                missed = (peptides[pos2] + peptides[pos2 + 1])
                if missed not in dict.keys():
                    dict[missed] = protein + "_" + str(pos_missed)
                else:
                    temp = str(dict[missed])
                    if protein + "_" + str(pos_missed) not in temp:
                        dict[missed] = temp + ";" + protein + "_" + str(pos_missed)

        # without N-terminal 'M'
        if seq.startswith('M'):
            pep_temp2 = re.sub(r'(?<=[RK])(?=[^P])', '\n', (seq[1:]), 0, re.DOTALL)
            pep2 = (pep_temp2.split())
            protein_length2 = 1
            pos3 = []
            peptides2 = []
            positions2 = []
            for p in pep2:
                start2 = protein_length2
                protein_length2 += len(p)
                pos3.append(start2)
            for i in zip(pep2, pos3):
                peptide2 = (i[0])
                peptides2.append(peptide2)
                position2 = (i[1])
                positions2.append(position2)
                if len(peptide2) >= 5:
                    if peptide2 not in dict.keys():
                        dict[peptide2] = protein + "_" + str(position2 + 1)
                    else:
                        temp = str(dict[peptide2])
                        if protein + "_" + str(position2+1) not in temp:
                            dict[peptide2] = temp + ";" + protein + "_" + str(position2 +1)


            count2 = 0
            for i in peptides2:
                pos4 = count2
                count2 += 1
                if pos4 <= (len(peptides2) - 5):
                    missed7 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3] +
                               peptides2[pos4 + 4])
                    missed6 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3])
                    missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                    missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                    pos_missed2 = (positions2[pos4])
                    if missed7 not in dict.keys():
                        dict[missed7] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed7])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed7] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed6 not in dict.keys():
                        dict[missed6] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed6])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed6] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed5 not in dict.keys():
                        dict[missed5] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed5])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed5] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed4 not in dict.keys():
                        dict[missed4] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed4])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed4] = temp + ";" + protein + "_" + str(pos_missed2 + 1)

                if pos4 <= (len(peptides2) - 4):
                    missed6 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3])
                    missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                    missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                    pos_missed2 = (positions2[pos4])
                    if missed6 not in dict.keys():
                        dict[missed6] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed6])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed6] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed5 not in dict.keys():
                        dict[missed5] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed5])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed5] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed4 not in dict.keys():
                        dict[missed4] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed4])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed4] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                if pos4 <= (len(peptides2) - 3):
                    missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                    missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                    pos_missed2 = (positions2[pos4])
                    if missed5 not in dict.keys():
                        dict[missed5] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed5])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed5] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                    if missed4 not in dict.keys():
                        dict[missed4] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed4])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed4] = temp + ";" + protein + "_" + str(pos_missed2 + 1)
                if pos4 <= (len(peptides2) - 2):
                    missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                    pos_missed2 = (positions2[pos4])
                    if missed4 not in dict.keys():
                        dict[missed4] = protein + "_" + str(pos_missed2 + 1)
                    else:
                        temp = str(dict[missed4])
                        if protein + "_" + str(pos_missed2 + 1) not in temp:
                            dict[missed4] = temp + ";" + protein + "_" + str(pos_missed2 + 1)

    return (dict)

RAP_DB = "database/fasta/RAP-DB.fasta"
RAP_DB_seq_dict = extract_peptides(RAP_DB)

MSU = "database/fasta/MSU.fasta"
MSU_seq_dict = extract_peptides(MSU)

UP = "database/fasta/Uniprot.fasta"
UP_seq_dict = extract_peptides(UP)

flr_filter=0.05

PXD000923 = "PXD000923/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df = pd.read_csv(PXD000923)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD000923_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD000923_peptides=df['Peptide_mod_pos'].unique()
PXD002222 = "PXD002222/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD002222)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD002222_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD002222_peptides=df['Peptide_mod_pos'].unique()
PXD002756 = "PXD002756/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD002756)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD002756_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD002756_peptides=df['Peptide_mod_pos'].unique()
PXD004705 = "PXD004705/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD004705)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD004705_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD004705_peptides=df['Peptide_mod_pos'].unique()
PXD004939 = "PXD004939/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD004939)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD004939_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD004939_peptides=df['Peptide_mod_pos'].unique()
PXD005241 = "PXD005241/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD005241)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD005241_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD005241_peptides=df['Peptide_mod_pos'].unique()
PXD012764 = "PXD012764/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD012764)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD012764_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD012764_peptides=df['Peptide_mod_pos'].unique()
PXD019291 = "PXD019291/FDR_0.01/binomial_peptidoform_collapsed_FLR.csv"
df=pd.read_csv(PXD019291)
df=df.loc[df['pA_q_value_BA']<=flr_filter]
df=df.sort_values(['Peptide_mod_pos','Binomial_final_score'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_mod_pos'),keep="last",inplace=False)
df['Score_FLR']=df['Binomial_final_score'].astype(str)+";"+df['pA_q_value_BA'].astype(str)+";"+df['USI']
PXD019291_dict=dict(zip(df.Peptide_mod_pos, df.Score_FLR))
PXD019291_peptides=df['Peptide_mod_pos'].unique()

PXD000923_peptides=list(np.array(PXD000923_peptides))
PXD002222_peptides=list(np.array(PXD002222_peptides))
PXD002756_peptides=list(np.array(PXD002756_peptides))
PXD004705_peptides=list(np.array(PXD004705_peptides))
PXD004939_peptides=list(np.array(PXD004939_peptides))
PXD005241_peptides=list(np.array(PXD005241_peptides))
PXD012764_peptides=list(np.array(PXD012764_peptides))
PXD019291_peptides=list(np.array(PXD019291_peptides))
peptide_list=PXD000923_peptides+PXD002222_peptides+PXD002756_peptides+PXD004705_peptides+PXD004939_peptides+PXD005241_peptides+PXD012764_peptides+PXD019291_peptides#+PXD001168_peptides
peptide_list=list(set(peptide_list))

peptide_list_no_mods = []
RAP_DB_list = []
MSU_list = []
rep_list =[]
UP_list = []
PXD000923_list=[]
PXD002222_list=[]
PXD002756_list=[]
PXD004705_list=[]
PXD004939_list=[]
PXD005241_list=[]
PXD012764_list=[]
PXD019291_list=[]
PXD000923_score_list=[]
PXD002222_score_list=[]
PXD002756_score_list=[]
PXD004705_score_list=[]
PXD004939_score_list=[]
PXD005241_score_list=[]
PXD012764_score_list=[]
PXD019291_score_list=[]
PXD000923_FLR_list=[]
PXD002222_FLR_list=[]
PXD002756_FLR_list=[]
PXD004705_FLR_list=[]
PXD004939_FLR_list=[]
PXD005241_FLR_list=[]
PXD012764_FLR_list=[]
PXD019291_FLR_list=[]
PXD000923_USI_list=[]
PXD002222_USI_list=[]
PXD002756_USI_list=[]
PXD004705_USI_list=[]
PXD004939_USI_list=[]
PXD005241_USI_list=[]
PXD012764_USI_list=[]
PXD019291_USI_list=[]


transcripts = pd.read_csv("SNP/all.locus_brief_info.7.0", sep="\t") #tsv file
transcripts = transcripts.loc[transcripts['is_representative']=="Y"]
transcripts_dict = dict(zip(transcripts.locus, transcripts.model))

SNP_file_loc="SNP/"

for i in peptide_list:
    peptidoform=i.rsplit("-",1)[0]
    peptide=re.sub(r'\[.*?\]','',peptidoform)
    peptide_list_no_mods.append(peptide)
    if peptide in RAP_DB_seq_dict.keys():
        temp_list=RAP_DB_seq_dict[peptide].split(";")
        final=""
        for temp in temp_list:
            pos=int(i.rsplit("-",1)[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"
        RAP_DB_list.append(final[:-1])
    else:
        RAP_DB_list.append("N/A")
    if peptide in MSU_seq_dict.keys():
        temp_list=MSU_seq_dict[peptide].split(";")
        final=""
        final_rep=""
        for temp in temp_list:
            #print(temp)
            pos=int(i.rsplit("-",1)[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"

            ##find representative - look up position in reference fasta
            seq=""
            try:
                rep=transcripts_dict[temp.split(".")[0]]
            except:
                rep=""
            if os.path.exists(SNP_file_loc+"output_SAAVs_June2023/"+rep+"_proteins_in_varieties.fasta.gz")==True:
                fasta=SNP_file_loc+"output_SAAVs_June2023/"+rep+"_proteins_in_varieties.fasta.gz"
            elif os.path.exists(SNP_file_loc + "output_SAAVs/output_SAAVs~/output_SAAVs/" + rep + "_proteins_in_varieties.fasta.gz") == True:
                fasta = SNP_file_loc + "output_SAAVs/output_SAAVs~/output_SAAVs/" + rep + "_proteins_in_varieties.fasta.gz"
            elif os.path.exists(SNP_file_loc + "out_saaps_all_canonical_MSU/" + rep + "_proteins_in_varieties.fasta.gz") == True:
                fasta = SNP_file_loc + "out_saaps_all_canonical_MSU/" + rep + "_proteins_in_varieties.fasta.gz"
            elif os.path.exists(SNP_file_loc + "out_saaps_msu_non_canonical/" + rep + "_proteins_in_varieties.fasta.gz") == True:
                fasta = SNP_file_loc + "out_saaps_msu_non_canonical/" + rep + "_proteins_in_varieties.fasta.gz"
            else:
                fasta="NOT FOUND"
                sub=""
                new_start = ""
                new_pos = ""
                #final_rep += "NA;"
                final_rep+=temp+"_"+str(pos)+";"
            if fasta!="NOT FOUND":
                with gzip.open(fasta, "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        if "reference" in record.id:
                            seq=record.seq
                new_start=seq.find(peptide)
                new_start=int(new_start)+1
                new_pos=new_start+int(i.rsplit("-",1)[1])-1
                final_rep+=rep+"_"+str(new_start)+"_"+str(new_pos)+";"


        MSU_list.append(final[:-1])
        rep_list.append(final_rep[:-1])

    else:
        MSU_list.append("N/A")
        rep_list.append("N/A")
    if peptide in UP_seq_dict.keys():
        temp_list=UP_seq_dict[peptide].split(";")
        final=""
        for temp in temp_list:
            pos=int(i.rsplit("-",1)[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"
        UP_list.append(final[:-1])
    else:
        UP_list.append("N/A")

    if i in PXD000923_peptides:
        PXD000923_list.append(1)
        PXD000923_score_list.append(PXD000923_dict[i].split(";")[0])
        PXD000923_FLR_list.append(PXD000923_dict[i].split(";")[1])
        PXD000923_USI_list.append(PXD000923_dict[i].split(";")[2])
    else:
        PXD000923_list.append(0)
        PXD000923_score_list.append("N/A")
        PXD000923_FLR_list.append("N/A")
        PXD000923_USI_list.append("N/A")
    if i in PXD002222_peptides:
        PXD002222_list.append(1)
        PXD002222_score_list.append(PXD002222_dict[i].split(";")[0])
        PXD002222_FLR_list.append(PXD002222_dict[i].split(";")[1])
        PXD002222_USI_list.append(PXD002222_dict[i].split(";")[2])
    else:
        PXD002222_list.append(0)
        PXD002222_score_list.append("N/A")
        PXD002222_FLR_list.append("N/A")
        PXD002222_USI_list.append("N/A")
    if i in PXD002756_peptides:
        PXD002756_list.append(1)
        PXD002756_score_list.append(PXD002756_dict[i].split(";")[0])
        PXD002756_FLR_list.append(PXD002756_dict[i].split(";")[1])
        PXD002756_USI_list.append(PXD002756_dict[i].split(";")[2])
    else:
        PXD002756_list.append(0)
        PXD002756_score_list.append("N/A")
        PXD002756_FLR_list.append("N/A")
        PXD002756_USI_list.append("N/A")
    if i in PXD004705_peptides:
        PXD004705_list.append(1)
        PXD004705_score_list.append(PXD004705_dict[i].split(";")[0])
        PXD004705_FLR_list.append(PXD004705_dict[i].split(";")[1])
        PXD004705_USI_list.append(PXD004705_dict[i].split(";")[2])
    else:
        PXD004705_list.append(0)
        PXD004705_score_list.append("N/A")
        PXD004705_FLR_list.append("N/A")
        PXD004705_USI_list.append("N/A")
    if i in PXD004939_peptides:
        PXD004939_list.append(1)
        PXD004939_score_list.append(PXD004939_dict[i].split(";")[0])
        PXD004939_FLR_list.append(PXD004939_dict[i].split(";")[1])
        PXD004939_USI_list.append(PXD004939_dict[i].split(";")[2])
    else:
        PXD004939_list.append(0)
        PXD004939_score_list.append("N/A")
        PXD004939_FLR_list.append("N/A")
        PXD004939_USI_list.append("N/A")
    if i in PXD005241_peptides:
        PXD005241_list.append(1)
        PXD005241_score_list.append(PXD005241_dict[i].split(";")[0])
        PXD005241_FLR_list.append(PXD005241_dict[i].split(";")[1])
        PXD005241_USI_list.append(PXD005241_dict[i].split(";")[2])
    else:
        PXD005241_list.append(0)
        PXD005241_score_list.append("N/A")
        PXD005241_FLR_list.append("N/A")
        PXD005241_USI_list.append("N/A")
    if i in PXD012764_peptides:
        PXD012764_list.append(1)
        PXD012764_score_list.append(PXD012764_dict[i].split(";")[0])
        PXD012764_FLR_list.append(PXD012764_dict[i].split(";")[1])
        PXD012764_USI_list.append(PXD012764_dict[i].split(";")[2])
    else:
        PXD012764_list.append(0)
        PXD012764_score_list.append("N/A")
        PXD012764_FLR_list.append("N/A")
        PXD012764_USI_list.append("N/A")
    if i in PXD019291_peptides:
        PXD019291_list.append(1)
        PXD019291_score_list.append(PXD019291_dict[i].split(";")[0])
        PXD019291_FLR_list.append(PXD019291_dict[i].split(";")[1])
        PXD019291_USI_list.append(PXD019291_dict[i].split(";")[2])
    else:
        PXD019291_list.append(0)
        PXD019291_score_list.append("N/A")
        PXD019291_FLR_list.append("N/A")
        PXD019291_USI_list.append("N/A")

df_final = pd.DataFrame()
df_final['Peptide_mod']=peptide_list
df_final['Peptide']=peptide_list_no_mods
df_final['RAP_DB']=RAP_DB_list
df_final['MSU']=MSU_list
df_final['UP']=UP_list
df_final['Representative_protein']=rep_list
df_final['PXD000923']=PXD000923_list
df_final['PXD002222']=PXD002222_list
df_final['PXD002756']=PXD002756_list
df_final['PXD004705']=PXD004705_list
df_final['PXD004939']=PXD004939_list
df_final['PXD005241']=PXD005241_list
df_final['PXD012764']=PXD012764_list
df_final['PXD019291']=PXD019291_list
df_final['PXD000923_score']=PXD000923_score_list
df_final['PXD002222_score']=PXD002222_score_list
df_final['PXD002756_score']=PXD002756_score_list
df_final['PXD004705_score']=PXD004705_score_list
df_final['PXD004939_score']=PXD004939_score_list
df_final['PXD005241_score']=PXD005241_score_list
df_final['PXD012764_score']=PXD012764_score_list
df_final['PXD019291_score']=PXD019291_score_list
df_final['PXD000923_FLR']=PXD000923_FLR_list
df_final['PXD002222_FLR']=PXD002222_FLR_list
df_final['PXD002756_FLR']=PXD002756_FLR_list
df_final['PXD004705_FLR']=PXD004705_FLR_list
df_final['PXD004939_FLR']=PXD004939_FLR_list
df_final['PXD005241_FLR']=PXD005241_FLR_list
df_final['PXD012764_FLR']=PXD012764_FLR_list
df_final['PXD019291_FLR']=PXD019291_FLR_list
df_final['PXD000923_USI']=PXD000923_USI_list
df_final['PXD002222_USI']=PXD002222_USI_list
df_final['PXD002756_USI']=PXD002756_USI_list
df_final['PXD004705_USI']=PXD004705_USI_list
df_final['PXD004939_USI']=PXD004939_USI_list
df_final['PXD005241_USI']=PXD005241_USI_list
df_final['PXD012764_USI']=PXD012764_USI_list
df_final['PXD019291_USI']=PXD019291_USI_list

for i in range(len(df_final)):
    flr1_count=0
    for database in ["PXD000923_FLR","PXD002222_FLR","PXD002756_FLR","PXD004705_FLR","PXD004939_FLR","PXD005241_FLR","PXD012764_FLR","PXD019291_FLR"]:
        if df_final.loc[i,database]!="N/A":
            if float(df_final.loc[i,database])<=0.01:
                flr1_count+=1
    if flr1_count>1:
        df_final.loc[i,'PTM_FLR_category']="Gold"
    elif flr1_count==1:
        df_final.loc[i,'PTM_FLR_category']="Silver"
    else:
        df_final.loc[i,'PTM_FLR_category']="Bronze"


USI=[]
df_temp = df_final.fillna(0)
df_temp= df_temp.replace("N/A",0)
for i in range(len(df_temp)):
    USI_dict = {"PXD000923": df_temp.loc[i, 'PXD000923_USI'],
                "PXD002222": df_temp.loc[i, 'PXD002222_USI'],
                "PXD002756": df_temp.loc[i, 'PXD002756_USI'],
                "PXD004705": df_temp.loc[i, 'PXD004705_USI'],
                "PXD004939": df_temp.loc[i, 'PXD004939_USI'],
                "PXD005241": df_temp.loc[i, 'PXD005241_USI'],
                "PXD012764": df_temp.loc[i, 'PXD012764_USI'],
                "PXD019291": df_temp.loc[i, 'PXD019291_USI']
                }

    score_dict = {float(df_temp.loc[i, 'PXD000923_score']): "PXD000923",
                  float(df_temp.loc[i, 'PXD002222_score']): "PXD002222",
                  float(df_temp.loc[i, 'PXD002756_score']): "PXD002756",
                  float(df_temp.loc[i, 'PXD004705_score']): "PXD004705",
                  float(df_temp.loc[i, 'PXD004939_score']): "PXD004939",
                  float(df_temp.loc[i, 'PXD005241_score']): "PXD005241",
                  float(df_temp.loc[i, 'PXD012764_score']): "PXD012764",
                  float(df_temp.loc[i, 'PXD019291_score']): "PXD019291",
                  }
    max_score = max(score_dict)
    PXD = score_dict[max_score]
    max_USI = USI_dict[PXD]
    USI.append(max_USI)
df_final["USI"]=USI

df_final=df_final.drop(columns=["PXD000923_USI","PXD002222_USI","PXD002756_USI","PXD004705_USI","PXD004939_USI","PXD005241_USI","PXD012764_USI","PXD019291_USI"])

df_final.to_csv("Rice_phosphosite_matrix_binomial_peptidoform_"+str(flr_filter)+"_w_protein-pos_scores_FLR.csv", index=False)

