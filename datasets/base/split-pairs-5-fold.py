# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:25:31 2023

@author: yangtingpeng
"""

import random 
seeds = [1234]

file="./pairs.txt"
with open(file,'r') as f:
    lines=f.readlines()

def save_file(file,data):
    with open(file,'w') as f:
        for line in data:
            f.write(line)

for seed in seeds:
    random.seed(seed)
    for i in range(1):
        random.shuffle(lines)
    indexs=list(range(len(lines)))
    folds=5
    step=len(lines)//folds
    for fold in range(folds):
        test=list(range(fold*step,(fold+1)*step))
        train=list(set(indexs)-set(test))

        train=[lines[i] for i in train]
        test=[lines[i] for i in test]

        save_file('./'+f'train_fold_{fold}.txt',train)
        save_file('./'+f'test_fold_{fold}.txt',test)

