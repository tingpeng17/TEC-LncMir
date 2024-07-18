#!/bin/bash
# $1: the lncRNA-miRNA pairs, $2 lncRNA fasta file, $3 miRNA fasta file
chmod 777 ./miranda
while IFS=$'\t' read -r lncrna mirna probability; do
    echo "lncRNA: $lncrna"
    echo "miRNA: $mirna"
    echo "--------------------"
    python prepare_for_miRanda.py --lncrna-fasta $2 --mirna-fasta $3 --lncrna $lncrna --mirna $mirna
    ./miranda mirna.fa lncrna.fa -en -1 -out $lncrna-$mirna-secondary-structure.txt
done < $1
rm lncrna.fa mirna.fa
