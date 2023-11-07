from __future__ import print_function, division
from Bio import SeqIO
import torch
import torch.utils.data
from datetime import datetime
import numpy as np

def log(m, file='None', timestamped=True, print_also=True):
    curr_time = f"[{datetime.now().strftime('%Y-%m-%d-%H:%M:%S')}] "
    log_string = f"{curr_time if timestamped else ''}{m}"
    if file =='None':
        print(log_string)
    else:
        with open(file, "a") as f:
            f.write(log_string+'\n')
        if print_also:
            print(log_string)

class PairedDataset(torch.utils.data.Dataset):
    """
    Dataset to be used by the PyTorch data loader for pairs of sequences and their labels.

    :param X0: List of first item in the pair
    :param X1: List of second item in the pair
    :param Y: List of labels
    """

    def __init__(self, X0, X1, Y, lnc_seq, mi_seq):
        self.X0 = X0
        self.X1 = X1
        self.Y = Y
        
        seqs_dict = {}
        seqs = list(SeqIO.parse(lnc_seq,format="fasta"))
        for x in seqs:
            id = str(x.id).split("|")[0]
            seq = str(x.seq)
            seqs_dict[id] = seq.upper()
        self.lnc_seqs_dict = seqs_dict
        
        seqs_dict = {}
        seqs = list(SeqIO.parse(mi_seq,format="fasta"))
        for x in seqs:
            id = str(x.id).split("|")[0]
            seq = str(x.seq)
            seqs_dict[id] = seq.upper()
        self.mi_seqs_dict = seqs_dict

    def __len__(self):
        return len(self.X0)

    def __getitem__(self, i):
        return self.lnc_seqs_dict[self.X0[i]], self.mi_seqs_dict[self.X1[i]], self.Y[i]


def collate_paired_sequences(args):
    """
    Collate function for PyTorch data loader.
    """
    x0 = [a[0] for a in args]
    x1 = [a[1] for a in args]
    y = [a[2] for a in args]
    return x0, x1, torch.stack(y, 0)
                    

def get_tokens(rna_list,base_number_dict):
    # rna_list: [ [">AGCT-AA","AGCUAA"], [">AGCT-AA","AGCUAA"] ]                
    for i in range(len(rna_list)):
        rna_list[i][1]=list(rna_list[i][1])
        for j in range(len(rna_list[i][1])):
            rna_list[i][1][j]=base_number_dict[rna_list[i][1][j]]
    return rna_list


def get_tokens_word(rna_list,base_number_dict):
    # rna_list: [ [">AGCT-AA","AGCUAA"], [">AGCT-AA","AGCUAA"] ]    
    one_word = int(np.log(len(base_number_dict))/np.log(4))
    for i in range(len(rna_list)):
        temp=[]
        for j in range(len(rna_list[i][1])//one_word):
            s=rna_list[i][1][(j*one_word):((j+1)*one_word)]
            temp.append(base_number_dict[s])
        rna_list[i][1]=temp
    return rna_list
