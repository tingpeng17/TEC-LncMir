# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 01:42:46 2023

@author: Tingpeng Yang

Evaluate a trained model.
"""
from __future__ import annotations
from Bio import SeqIO
import argparse
import datetime
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score,precision_recall_curve,roc_auc_score,roc_curve
from sklearn.metrics import confusion_matrix
from tqdm import tqdm
from utils import log,get_tokens,get_tokens_word
matplotlib.use("Agg")
def add_args(parser):
    # Create parser for command line utility.
    parser.add_argument("--model", help="Trained prediction model", required=True)
    parser.add_argument("--test-pair", required=True, help="the url of the test pairs file", type=str)
    parser.add_argument("--test-lnc-seq", required=True, help="the url of the lncrna fasta file for lncrnas in test pairs", type=str)
    parser.add_argument("--test-mi-seq", required=True, help="the url of the mirna fasta file for mirnas in test pairs", type=str)
    parser.add_argument("-o", "--outdir", help="Output file to write results")
    parser.add_argument("-d", "--device", type=int, default=0, help="Compute device to use")
    return parser

def plot_eval_predictions(labels, predictions,output, path="figure"):
    """
    Plot histogram of positive and negative predictions, precision-recall curve, and receiver operating characteristic curve.
    :param y: Labels
    :type y: np.ndarray
    :param phat: Predicted probabilities
    :type phat: np.ndarray
    :param path: File prefix for plots to be saved to [default: figure]
    :type path: str
    """
    pos_phat = predictions[labels == 1]
    neg_phat = predictions[labels == 0]
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle("Distribution of Predictions")
    ax1.hist(pos_phat)
    ax1.set_xlim(0, 1)
    ax1.set_title("Positive")
    ax1.set_xlabel("p-hat")
    ax2.hist(neg_phat)
    ax2.set_xlim(0, 1)
    ax2.set_title("Negative")
    ax2.set_xlabel("p-hat")
    plt.savefig(path + "phat_dist.svg")
    plt.close()
    precision, recall, pr_thresh = precision_recall_curve(labels, predictions)
    aupr = average_precision_score(labels, predictions)
    log(f"AUPR: {aupr}",file=output)
    plt.step(recall, precision, color="b", alpha=0.2, where="post")
    plt.fill_between(recall, precision, step="post", alpha=0.2, color="b")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title("Precision-Recall (AUPR: {:.3})".format(aupr))
    plt.savefig(path + "aupr.svg")
    plt.close()
    fpr, tpr, roc_thresh = roc_curve(labels, predictions)
    auroc = roc_auc_score(labels, predictions)
    log(f"AUROC: {auroc}",file=output)
    plt.step(fpr, tpr, color="b", alpha=0.2, where="post")
    plt.fill_between(fpr, tpr, step="post", alpha=0.2, color="b")
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title("Receiver Operating Characteristic (AUROC: {:.3})".format(auroc))
    plt.savefig(path + "auroc.svg")
    plt.close()

def main(args):
    """
    Run model evaluation from arguments.
    """
    if args.outdir is None:
        outPath = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")+'-'
    else:
        outPath = args.outdir
    outFile = open(outPath + "predictions.tsv", "w+")
    output=outPath+'log.txt'
    # Set Device
    device = args.device
    output=args.outdir+'log.txt' # the log_file
    if device == -1 :
        args.device='cpu'
        device = args.device
        log(
            "Using cpu",
            file=output,
            print_also=True,
        )
    elif device > -1 :
        log(
            f"Using CUDA device {device} - {torch.cuda.get_device_name(device)}",
            file=output,
            print_also=True,
        )
    else:
        log(
            "device must be an integer >= -1!",
            file=output,
            print_also=True,
        )
        sys.exit('device must be an integer >= -1!')
    # Load Model
    model_path = args.model
    model = torch.load(model_path).to(args.device)
    model.use_cuda = True
    #for name,parameters in model.named_parameters():
    #    print(name,':',parameters.size())
    # Load Pairs
    test_fi = args.test_pair
    test_df = pd.read_csv(test_fi, sep="\t", header=None)
    
    test_lnc_seqs_dict = {}
    seqs = list(SeqIO.parse(args.test_lnc_seq,format="fasta"))
    for x in seqs:
        id = str(x.id).split("|")[0]
        seq = str(x.seq).upper()
        test_lnc_seqs_dict[id] = seq
        
    test_mi_seqs_dict = {}
    seqs = list(SeqIO.parse(args.test_mi_seq,format="fasta"))
    for x in seqs:
        id = str(x.id).split("|")[0]
        seq = str(x.seq).upper()
        test_mi_seqs_dict[id] = seq
    
    
    
    model.eval()
    base_number_dict_mirna={"A":1,"G":2,"C":3,"U":4}
    base_number_dict_lnc={}
    k=0
    if model.one_word ==1 :
        base_number_dict_lnc={"A":1,"G":2,"C":3,"U":4}
    elif model.one_word ==2 :
        for x in "AUCG":
            for y in "AUCG":
                k=k+1
                s=x+y
                base_number_dict_lnc[s]=k
    elif model.one_word ==3 :
        for x in "AUCG":
            for y in "AUCG":
                for z in "AUCG":
                    k=k+1
                    s=x+y+z
                    base_number_dict_lnc[s]=k
    elif model.one_word ==4 :
        for x in "AUCG":
            for y in "AUCG":
                for z in "AUCG":
                    for p in "AUCG":
                        k=k+1
                        s=x+y+z+p
                        base_number_dict_lnc[s]=k
    elif model.one_word ==5 :
        for x in "AUCG":
            for y in "AUCG":
                for z in "AUCG":
                    for p in "AUCG":
                        for q in "AUCG":
                            k=k+1
                            s=x+y+z+p+q
                            base_number_dict_lnc[s]=k
    elif model.one_word ==6 :
        for x in "AUCG":
            for y in "AUCG":
                for z in "AUCG":
                    for p in "AUCG":
                        for q in "AUCG":
                            for m in "AUCG":
                                k=k+1
                                s=x+y+z+p+q+m
                                base_number_dict_lnc[s]=k
    elif model.one_word ==7 :
        for x in "AUCG":
            for y in "AUCG":
                for z in "AUCG":
                    for p in "AUCG":
                        for q in "AUCG":
                            for m in "AUCG":
                                for n in "AUCG":
                                    k=k+1
                                    s=x+y+z+p+q+m+n
                                    base_number_dict_lnc[s]=k       
    with torch.no_grad():
        phats = []
        labels = []
        for _, (n0, n1, label) in tqdm(test_df.iterrows(), total=len(test_df), desc="Predicting pairs"):
            lncrnas = list(set([test_lnc_seqs_dict[n0]]))
            lncrnas=[[i,i.replace('-','').replace('>','').replace("T", "U")] for i in lncrnas]
            mirnas = list(set([test_mi_seqs_dict[n1]]))
            mirnas=[[i,i.replace('-','').replace('>','').replace("T", "U")] for i in mirnas]
            rna_list=get_tokens(mirnas,base_number_dict_mirna)+get_tokens_word(lncrnas,base_number_dict_lnc)
            for i in range(len(rna_list)):
                rna_list[i][1]=torch.LongTensor(rna_list[i][1]).to(args.device)
            embeddings={rna_list[i][0]:rna_list[i][1] for i in range(len(rna_list))}
            p0 = embeddings[test_lnc_seqs_dict[n0]].reshape(1,-1,1)
            p1 = embeddings[test_mi_seqs_dict[n1]].reshape(1,-1,1)
            p0=p0.to(args.device)
            p1=p1.to(args.device)
            pred = model.map_predict(p0, p1)[1].item()
            phats.append(pred)
            labels.append(label)
            outFile.write(f"{n0}\t{n1}\t{label}\t{pred:.5}\n")
    if len(set(labels))==1:
        phats = np.array(phats)
        labels = np.array(labels)
        predicted_labels=phats 
        predicted_labels[predicted_labels>=0.5]=1
        predicted_labels[predicted_labels<0.5]=0
        accuracy=accuracy_score(labels, predicted_labels)
        log("accuracy: {:.4f}".format(accuracy),file=output)
        outFile.close()
    else:
        phats = np.array(phats)
        labels = np.array(labels)
        plot_eval_predictions(labels, phats, output,outPath)
        predicted_labels=phats 
        predicted_labels[predicted_labels>=0.5]=1
        predicted_labels[predicted_labels<0.5]=0
        tn, fp, fn, tp = confusion_matrix(labels, predicted_labels).ravel()
        accuracy = (tp+tn)/(tp+tn+fp+fn)
        sensitivity = tp/(tp+fn) 
        specificity = tn/(tn+fp)
        ppv = tp/(tp+fp)
        npv = tn/(tn+fn)
        f1score = 2*tp/(2*tp+fp+fn)
        mcc=(tp*tn-fp*fn)/(((tp+fn)*(tp+fp)*(tn+fp)*(tn+fn))**0.5)
        log("accuracy: {:.4f},sensitivity: {:.4f},specificity: {:.4f},ppv: {:.4f},npv: {:.4f},F1 score: {:.4f},mcc: {:.4f}".format(accuracy,sensitivity,specificity,ppv,npv,f1score,mcc),file=output)
        outFile.close()
        
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser = add_args(parser)
    main(parser.parse_args())