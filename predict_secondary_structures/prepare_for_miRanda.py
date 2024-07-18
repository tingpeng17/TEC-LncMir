from Bio import SeqIO
import argparse

def add_args(parser):
    parser.add_argument("--mirna-fasta", required=True, help="the url of the mirna file", type=str)
    parser.add_argument("--lncrna-fasta", required=True, help="the url of the lncrna file", type=str)
    parser.add_argument("--lncrna", required=True, help="the name of lncrna", type=str)
    parser.add_argument("--mirna", required=True, help="the name of mirna", type=str)
    return parser

parser = argparse.ArgumentParser()
parser = add_args(parser)
args =parser.parse_args()


lncs_dict={}
seqs = list(SeqIO.parse(args.lncrna_fasta,format="fasta"))
for x in seqs:
    id = str(x.id).split("|")[0]
    seq = str(x.seq).upper().replace('T','U')
    lncs_dict[id] = seq
mirnas_dict={}
seqs = list(SeqIO.parse(args.mirna_fasta,format="fasta"))
for x in seqs:
    id = str(x.id).split("|")[0]
    seq = str(x.seq).upper().replace('T','U')
    mirnas_dict[id] = seq
    
with open('lncrna.fa','w') as f1:
    with open('mirna.fa','w') as f2:
        f1.write(f"{'>'+args.lncrna}\n{lncs_dict[args.lncrna]}\n")
        f2.write(f"{'>'+args.mirna}\n{mirnas_dict[args.mirna]}\n")