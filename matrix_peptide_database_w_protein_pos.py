import pandas as pd
from Bio import SeqIO
import re
import numpy as np

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

RAP_DB = "C:/Users/krams/Dropbox/PTMExchange/Rice/database/fasta/RAP-DB.fasta"
RAP_DB_seq_dict = extract_peptides(RAP_DB)

MSU = "C:/Users/krams/Dropbox/PTMExchange/Rice/database/fasta/MSU.fasta"
MSU_seq_dict = extract_peptides(MSU)

UP = "C:/Users/krams/Dropbox/PTMExchange/Rice/database/fasta/Uniprot.fasta"
UP_seq_dict = extract_peptides(UP)

flr_filter=0.05

PXD000923 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD000923/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df = pd.read_csv(PXD000923)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD000923_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD000923_peptides=df['Peptide_pos'].unique()
PXD002222 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD002222/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD002222)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD002222_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD002222_peptides=df['Peptide_pos'].unique()
PXD002756 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD002756/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD002756)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD002756_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD002756_peptides=df['Peptide_pos'].unique()
PXD004705 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD004705/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD004705)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD004705_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD004705_peptides=df['Peptide_pos'].unique()
PXD004939 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD004939/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD004939)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD004939_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD004939_peptides=df['Peptide_pos'].unique()
PXD005241 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD005241/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD005241)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD005241_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD005241_peptides=df['Peptide_pos'].unique()
PXD012764 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD012764/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD012764)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD012764_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD012764_peptides=df['Peptide_pos'].unique()
PXD019291 = "C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/PXD019291/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD019291)
df=df.loc[df['pAla_q_value']<=flr_filter]
df['Peptide_pos']=df['Peptide']+"-"+df['PTM positions'].astype(str)
df=df.sort_values(['Peptide_pos','PTM_final_prob'],ascending=[True,True])
df=df.drop_duplicates(subset=('Peptide_pos'),keep="last",inplace=False)
df['Score_FLR']=df['PTM_final_prob'].astype(str)+";"+df['pAla_q_value'].astype(str)
PXD019291_dict=dict(zip(df.Peptide_pos, df.Score_FLR))
PXD019291_peptides=df['Peptide_pos'].unique()

PXD000923_peptides=list(np.array(PXD000923_peptides))
PXD002222_peptides=list(np.array(PXD002222_peptides))
PXD002756_peptides=list(np.array(PXD002756_peptides))
PXD004705_peptides=list(np.array(PXD004705_peptides))
PXD004939_peptides=list(np.array(PXD004939_peptides))
PXD005241_peptides=list(np.array(PXD005241_peptides))
PXD012764_peptides=list(np.array(PXD012764_peptides))
PXD019291_peptides=list(np.array(PXD019291_peptides))
peptide_list=PXD000923_peptides+PXD002222_peptides+PXD002756_peptides+PXD004705_peptides+PXD004939_peptides+PXD005241_peptides+PXD012764_peptides+PXD019291_peptides
peptide_list=list(set(peptide_list))

RAP_DB_list = []
MSU_list = []
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

for i in peptide_list:
    if i.split("-")[0] in RAP_DB_seq_dict.keys():
        temp_list=RAP_DB_seq_dict[i.split("-")[0]].split(";")
        final=""
        for temp in temp_list:
            pos=int(i.split("-")[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"
        RAP_DB_list.append(final[:-1])
    else:
        RAP_DB_list.append("N/A")
    if i.split("-")[0] in MSU_seq_dict.keys():
        temp_list=MSU_seq_dict[i.split("-")[0]].split(";")
        final=""
        for temp in temp_list:
            pos=int(i.split("-")[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"
        MSU_list.append(final[:-1])
    else:
        MSU_list.append("N/A")
    if i.split("-")[0] in UP_seq_dict.keys():
        temp_list=UP_seq_dict[i.split("-")[0]].split(";")
        final=""
        for temp in temp_list:
            pos=int(i.split("-")[1])+int(temp.split("_")[-1])-1
            final+=temp+"_"+str(pos)+";"
        UP_list.append(final[:-1])
    else:
        UP_list.append("N/A")

    if i in PXD000923_peptides:
        PXD000923_list.append(1)
        PXD000923_score_list.append(PXD000923_dict[i].split(";")[0])
        PXD000923_FLR_list.append(PXD000923_dict[i].split(";")[1])
    else:
        PXD000923_list.append(0)
        PXD000923_score_list.append("N/A")
        PXD000923_FLR_list.append("N/A")
    if i in PXD002222_peptides:
        PXD002222_list.append(1)
        PXD002222_score_list.append(PXD002222_dict[i].split(";")[0])
        PXD002222_FLR_list.append(PXD002222_dict[i].split(";")[1])
    else:
        PXD002222_list.append(0)
        PXD002222_score_list.append("N/A")
        PXD002222_FLR_list.append("N/A")
    if i in PXD002756_peptides:
        PXD002756_list.append(1)
        PXD002756_score_list.append(PXD002756_dict[i].split(";")[0])
        PXD002756_FLR_list.append(PXD002756_dict[i].split(";")[1])
    else:
        PXD002756_list.append(0)
        PXD002756_score_list.append("N/A")
        PXD002756_FLR_list.append("N/A")
    if i in PXD004705_peptides:
        PXD004705_list.append(1)
        PXD004705_score_list.append(PXD004705_dict[i].split(";")[0])
        PXD004705_FLR_list.append(PXD004705_dict[i].split(";")[1])
    else:
        PXD004705_list.append(0)
        PXD004705_score_list.append("N/A")
        PXD004705_FLR_list.append("N/A")
    if i in PXD004939_peptides:
        PXD004939_list.append(1)
        PXD004939_score_list.append(PXD004939_dict[i].split(";")[0])
        PXD004939_FLR_list.append(PXD004939_dict[i].split(";")[1])
    else:
        PXD004939_list.append(0)
        PXD004939_score_list.append("N/A")
        PXD004939_FLR_list.append("N/A")
    if i in PXD005241_peptides:
        PXD005241_list.append(1)
        PXD005241_score_list.append(PXD005241_dict[i].split(";")[0])
        PXD005241_FLR_list.append(PXD005241_dict[i].split(";")[1])
    else:
        PXD005241_list.append(0)
        PXD005241_score_list.append("N/A")
        PXD005241_FLR_list.append("N/A")
    if i in PXD012764_peptides:
        PXD012764_list.append(1)
        PXD012764_score_list.append(PXD012764_dict[i].split(";")[0])
        PXD012764_FLR_list.append(PXD012764_dict[i].split(";")[1])
    else:
        PXD012764_list.append(0)
        PXD012764_score_list.append("N/A")
        PXD012764_FLR_list.append("N/A")
    if i in PXD019291_peptides:
        PXD019291_list.append(1)
        PXD019291_score_list.append(PXD019291_dict[i].split(";")[0])
        PXD019291_FLR_list.append(PXD019291_dict[i].split(";")[1])
    else:
        PXD019291_list.append(0)
        PXD019291_score_list.append("N/A")
        PXD019291_FLR_list.append("N/A")

df_final = pd.DataFrame()
df_final['Peptide']=peptide_list
df_final['RAP_DB']=RAP_DB_list
df_final['MSU']=MSU_list
df_final['UP']=UP_list
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

for i in range(len(df_final)):
    flr1_count=0
    for database in ["PXD000923_FLR","PXD002222_FLR","PXD002756_FLR","PXD004705_FLR","PXD004939_FLR","PXD005241_FLR","PXD012764_FLR","PXD019291_FLR"]:
        #df_final[database]=df_final[database].astype(float)
        #df_final=df_final.fillna("-")
        if df_final.loc[i,database]!="N/A":
            if float(df_final.loc[i,database])<=0.01:
                flr1_count+=1
    if flr1_count>1:
        df_final.loc[i,'PTM_FLR_category']="Gold"
    elif flr1_count==1:
        df_final.loc[i,'PTM_FLR_category']="Silver"
    else:
        df_final.loc[i,'PTM_FLR_category']="Bronze"

#df_final = pd.DataFrame(list(zip(peptide_list,RAP_DB_list,MSU_list,UP_list,PXD000923_list,PXD002222_list,PXD002756_list,PXD004705_list,PXD004939_list,PXD005241_list,
                                 #PXD012764_list,PXD019291_list)), columns = ['Peptide', 'RAP_DB', 'MSU', 'UP',"PXD000923","PXD002222","PXD002756","PXD004705","PXD004939","PXD005241","PXD012764","PXD019291"])
df_final.to_csv("C:/Users/krams/Dropbox/PTMExchange/Rice/New_build_ID/Rice_phosphosite_matrix_updated"+str(flr_filter)+"_w_protein-pos_scores_FLR.csv", index=False)
