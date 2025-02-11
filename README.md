# TEC_LncMir
TEC_LncMir is a deep learning model to predict interactions between lncRNAs and microRNAs based on the Transformer Encoder and convolutional neural networks (CNNs). TEC_LncMir uses the sequences of lncRNAs and microRNAs as the inputs and predicts the interaction probabilities between lncRNAs and microRNAs. If you use TEC_LncMir in your work, please cite the following publication:
Tingpeng Yang, Yonghong He, Yu Wang, Introducing TEC-LncMir for prediction of lncRNA-miRNA interactions through deep learning of RNA sequences, Briefings in Bioinformatics, Volume 26, Issue 1, January 2025, bbaf046, https://doi.org/10.1093/bib/bbaf046
# The usage of our code
## Preparation: 
Enter the TEC_LncMir folder
### Setup the enviroment

```
conda env create --file environment.yml
```

### Activate the enviroment

```
conda activate TEC_LncMir
```

### Prepare the folders for code running

```
mkdir output
mkdir output/model
mkdir output/evaluate
mkdir output/predict
```

### Prepare the datasets

```
cd datasets/base
python split-pairs-5-fold.py
python extend_pairs.py
cd ../large
unzip outLncRNA.zip
python split-pairs-5-fold.py
python extend_pairs.py
cd ../../
```

The following commands are some examples to help users easily use TEC_LncMir.

## Train a model from scratch:

```
python ./code/train.py --train-pair ./datasets/large/train_fold_0.txt --train-mi-seq ./datasets/large/homo_mature_mirna.fa --train-lnc-seq ./datasets/large/outLncRNA.fa --outdir ./output/model/ --valid-pair ./datasets/large/test_fold_0.txt --valid-mi-seq ./datasets/large/homo_mature_mirna.fa --valid-lnc-seq ./datasets/large/outLncRNA.fa
```

## Train a model from a checkpoint:

```
python ./code/train.py --train-pair ./datasets/large/train_fold_0.txt --train-mi-seq ./datasets/large/homo_mature_mirna.fa --train-lnc-seq ./datasets/large/outLncRNA.fa --outdir ./output/model/ --valid-pair ./datasets/large/test_fold_0.txt --valid-mi-seq ./datasets/large/homo_mature_mirna.fa --valid-lnc-seq ./datasets/large/outLncRNA.fa --checkpoint ./model_weights/large/control/fold0/model/best_model/
```

## Evaluate a pretrained model

```
python ./code/evaluate.py --test-pair ./datasets/large/test_fold_0.txt --test-mi-seq ./datasets/large/homo_mature_mirna.fa --test-lnc-seq ./datasets/large/outLncRNA.fa --device 0 --outdir ./model_weights/large/control/fold0/evaluate/ --model ./model_weights/large/control/fold0/model/best_model/model.sav
```

## Predict interactions between lncRNAs and miRNAs

Users need to prepare their data in the form of [pairs-need-predict.txt](pairs-need-predict.txt), and then run the following command.

```
python ./code/predict.py --pairs ./pairs-need-predict.txt --mi-seq ./datasets/large/homo_mature_mirna.fa --lnc-seq ./datasets/large/outLncRNA.fa --device 0 -o ./model_weights/large/control/fold0/predict/ --model ./model_weights/large/control/fold0/model/best_model/model.sav
```

The results will be shown in the dir ./model_weights/large/control/fold0/predict/.

The file predict.tsv shows all the predictions of TEC-LncMir.

The file predict-positive.tsv shows the positive interactions predicted by TEC-LncMir.

## Predict the secondary structures of the lncRNA-miRNA interactions

```
cd predict_secondary_structures
bash main.sh ../pairs-need-predict.txt ../datasets/large/outLncRNA.fa ../datasets/large/homo_mature_mirna.fa
```

The results will be shown in the dir predict_secondary_structures.

The secondary structures of the interactions for each pair of lncRNA and miRNA will be shown as {lncrna_name}-{mirna_name}-secondary-structure.txt
