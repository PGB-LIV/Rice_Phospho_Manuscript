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

PXD000923 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD000923/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df = pd.read_csv(PXD000923)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD000923_peptides=df['Peptide'].unique()
PXD002222 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD002222/6_files/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD002222)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD002222_peptides=df['Peptide'].unique()
PXD002756 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD002756/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD002756)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD002756_peptides=df['Peptide'].unique()
PXD004705 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD004705/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD004705)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD004705_peptides=df['Peptide'].unique()
PXD004939 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD004939/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD004939)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD004939_peptides=df['Peptide'].unique()
PXD005241 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD005241/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD005241)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD005241_peptides=df['Peptide'].unique()
PXD012764 = "C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/PXD012764/comet_tryptic/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD012764)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD012764_peptides=df['Peptide'].unique()
PXD019291 = "C:/Users/krams/Dropbox/PTMExchange/Rice/PXD019291/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_pAla.csv"
df=pd.read_csv(PXD019291)
df=df.loc[df['pAla_q_value']<=flr_filter]
PXD019291_peptides=df['Peptide'].unique()

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

for i in peptide_list:
    if i in RAP_DB_seq_dict.keys():
        RAP_DB_list.append(RAP_DB_seq_dict[i])
    else:
        RAP_DB_list.append("N/A")
    if i in MSU_seq_dict.keys():
        MSU_list.append(MSU_seq_dict[i])
    else:
        MSU_list.append("N/A")
    if i in UP_seq_dict.keys():
        UP_list.append(UP_seq_dict[i])
    else:
        UP_list.append("N/A")
    if i in PXD000923_peptides:
        PXD000923_list.append(1)
    else:
        PXD000923_list.append(0)
    if i in PXD002222_peptides:
        PXD002222_list.append(1)
    else:
        PXD002222_list.append(0)
    if i in PXD002756_peptides:
        PXD002756_list.append(1)
    else:
        PXD002756_list.append(0)
    if i in PXD004705_peptides:
        PXD004705_list.append(1)
    else:
        PXD004705_list.append(0)
    if i in PXD004939_peptides:
        PXD004939_list.append(1)
    else:
        PXD004939_list.append(0)
    if i in PXD005241_peptides:
        PXD005241_list.append(1)
    else:
        PXD005241_list.append(0)
    if i in PXD012764_peptides:
        PXD012764_list.append(1)
    else:
        PXD012764_list.append(0)
    if i in PXD019291_peptides:
        PXD019291_list.append(1)
    else:
        PXD019291_list.append(0)

df_final = pd.DataFrame(list(zip(peptide_list,RAP_DB_list,MSU_list,UP_list,PXD000923_list,PXD002222_list,PXD002756_list,PXD004705_list,PXD004939_list,PXD005241_list,
                                 PXD012764_list,PXD019291_list)), columns = ['Peptide', 'RAP_DB', 'MSU', 'UP',"PXD000923","PXD002222","PXD002756","PXD004705","PXD004939","PXD005241","PXD012764","PXD019291"])
df_final.to_csv("C:/Users/krams/Dropbox/PTMExchange/Rice/Eric/Rice_phosphosite_matrix_updated"+str(flr_filter)+".csv", index=False)
