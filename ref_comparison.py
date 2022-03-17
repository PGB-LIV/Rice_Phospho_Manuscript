# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:46:05 2022

@author: krams
"""


import os
import gzip
from Bio import SeqIO
import pandas as pd

SNP_file_loc="C:/Users/krams/Dropbox/PTMExchange/Rice/SNP/"

p="LOC_Os02g58110.1_442_444"
pos=466

record_all=[]
pos_all=[]

if os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
    fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
    file=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
elif os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
    fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
    file=NP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
else:
    print("Something is wrong")
if fasta!="NOT_FOUND":
    with gzip.open(fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if "reference" in record.id:
                seq=record.seq
else:
    seq="NA"
print(seq[pos-1],file)

if os.path.exists(SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
    fasta=SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
    file=SNP_file_loc+"out_saaps_msu_non_canonical/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
elif os.path.exists(SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz")==True:
    fasta=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
    file=SNP_file_loc+"out_saaps_all_canonical_MSU/"+p.rsplit("_",2)[0]+"_different_proteins_in_varieties.fasta.gz"
else:
    print("Something is wrong")
if fasta!="NOT_FOUND":
    with gzip.open(fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_all.append(record.id)
            pos_all.append(record.seq[pos-1])
            if "reference" in record.id:
                seq=record.seq
else:
    seq="NA"
print(seq[pos-1],file)

print(seq)

df=pd.DataFrame({"ID":record_all,pos:pos_all})
df.to_csv(SNP_file_loc+"seq_temp_"+p+"_pos_"+str(pos)+".csv")