# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:25:31 2023

@author: yangtingpeng
"""
kmer=4
lnc_fasta="./outLncRNA.fa"
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
lncrna = {}
list_lncrna = list(SeqIO.parse(lnc_fasta,
                               format="fasta"))
for x in list_lncrna:
    id = str(x.id).split("|")[0]
    seq = str(x.seq)
    lncrna[id] = seq

    for k in range(1,kmer):
        lncrna[id+f'-extend{k}'] = seq[k:]

records = []
for seq_id, sequence in lncrna.items():
    seq = Seq(sequence)
    record = SeqRecord(seq, id=seq_id)
    records.append(record)
SeqIO.write(records, lnc_fasta+'.extend', "fasta")

files=['train_fold_0.txt','train_fold_1.txt','train_fold_2.txt','train_fold_3.txt','train_fold_4.txt',
      #'test_fold_0.txt','test_fold_1.txt','test_fold_2.txt','test_fold_3.txt','test_fold_4.txt'
      ]
for file in files:
    with open(file,'r') as f:
        lines=f.readlines()
    lines=[line.strip().split('\t') for line in lines]
    outfile='extend_'+file
    re=[]
    re.append(lines)
    for k in range(1,kmer):
        temp=[[line[0]+f'-extend{k}',line[1],line[2]] for line in lines]
        re.append(temp)
    with open(outfile,'w') as f:
        for i in range(len(re[0])):
            for j in range(len(re)):
                f.write(f'{re[j][i][0]}\t{re[j][i][1]}\t{re[j][i][2]}\n')