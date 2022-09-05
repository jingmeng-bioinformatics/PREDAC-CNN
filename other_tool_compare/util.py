from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os, sys
import pandas as pd
import numpy as np
import scipy as sp
import random
import math
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt

from model import fasta_to_csv
from model import calculate_label
from model import generate_feature
from model import strain_selection
from model import replace_uncertain_amino_acids
from model import train_test_split_data

def distance_mutation(distance_input, seq_input):
    num_mut_list = []
    print("distance_mutation")
    for i in range(0, distance_input.shape[0]):
        print("---{:.2f}---".format(100*i/distance_input.shape[0]),end="\r")
        mutation_count = 0
        strain_1 = []
        strain_2 = []
        for j in range(0, seq_input.shape[0]):
            if seq_input['description'].iloc[j] == distance_input['Strain1'].iloc[i]:
                strain_1 = seq_input['seq'].iloc[j]
            if seq_input['description'].iloc[j] == distance_input['Strain2'].iloc[i]:
                strain_2 = seq_input['seq'].iloc[j]
        mutation_count = sum(0 if c1 == c2 else 1 for c1, c2 in zip(strain_1, strain_2))
        num_mut_list.append(mutation_count)
    return num_mut_list
        

def get_confusion_matrix(x_true, y_true, subtype):
    TP, FP, TN, FN = 0, 0, 0, 0
    threshold = 0
    #optimized threshold (with best performance)
    if subtype == 'H1N1':
        threshold = 11
    elif subtype == 'H3N2':
        threshold = 9
    elif subtype == 'H5N1':
        threshold = 12
        
    for i in range(len(y_true)):
        if x_true[i] <= threshold:  #<=
            if y_true[i] < 4:
                TN = TN + 1
            else:
                FP = FP + 1
        else:
            if y_true[i] >= 4:
                TP = TP + 1
            else:
                FN = FN + 1
                
    conf_matrix = [
        [FP, TP],
        [TN, FN]
    ]

    return conf_matrix

def get_confusion_matrix_test(x_true, y_true, subtype):
    threshold = 9
    ans_pre = []
    ans_valid = []
    #optimized threshold (with best performance)
    if subtype == 'H1N1':
        threshold = 11
    elif subtype == 'H3N2':
        threshold = 9
    elif subtype == 'H5N1':
        threshold = 12
    
    for i in range(len(x_true)):
        if(x_true[i] <=threshold):
            ans_pre.append(0)
        else:
            ans_pre.append(1)
        if(y_true[i] < 4):
            ans_valid.append(0)
        else:
            ans_valid.append(1)
    return ans_pre,ans_valid



def get_accuracy(conf_matrix):
    """
    Calculates accuracy metric from the given confusion matrix.
    """
    TP, FP, FN, TN = conf_matrix[0][1], conf_matrix[0][0], conf_matrix[1][1], conf_matrix[1][0]
    return (TP + TN) / (TP + FP + FN + TN)

def get_precision(conf_matrix):
    """
    Calculates precision metric from the given confusion matrix.
    """
    TP, FP = conf_matrix[0][1], conf_matrix[0][0]

    if TP + FP > 0:
        return TP / (TP + FP)
    else:
        return 0

def get_recall(conf_matrix):
    """
    Calculates recall metric from the given confusion matrix.
    """
    TP, FN = conf_matrix[0][1], conf_matrix[1][1]

    if TP + FN > 0:
        return TP / (TP + FN)
    else:
        return 0
    
def get_f1score(conf_matrix):
    """
    Calculates f1-score metric from the given confusion matrix.
    """
    p = get_precision(conf_matrix)
    r = get_recall(conf_matrix)

    if p + r > 0:
        return 2 * p * r / (p + r)
    else:
        return 0
    
def get_mcc(conf_matrix):
    """
    Calculates Matthew's Correlation Coefficient metric from the given confusion matrix.
    """
    TP, FP, FN, TN = conf_matrix[0][1], conf_matrix[0][0], conf_matrix[1][1], conf_matrix[1][0]
    if TP + FP > 0 and TP + FN > 0 and TN + FP > 0 and TN + FN > 0:
        return (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    else:
        return 0



































