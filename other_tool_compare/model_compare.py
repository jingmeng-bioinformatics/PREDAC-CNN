import os, sys
import numpy as np
import pandas as pd
import torch
import warnings
sys.path.append(os.path.abspath("/Share/home/goldwaterkk/AntigenPre/tool-compare/IAV-CNN-years/code_years"))
#sys.path.append(os.path.abspath("/content/drive/My Drive/Colab Notebooks/bioinformatics/code"))
from model import *
from train_cnn import train_cnn,test_cnn,test_cnn_ROC
from util import get_confusion_matrix
from util import get_confusion_matrix_test
from util import distance_mutation
from util import get_accuracy
from util import get_precision
from util import get_recall
from util import get_f1score
from util import get_mcc
warnings.filterwarnings('ignore')
def H3N2_Npy(years1,years2,types):
    years_data = []
    for i in range(years1,years2+1):
        try:
            years_data.append(np.load("./training_data/IVA-CNN/{}_{}.npy".format(str(i),types)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data),axis=0)
    #np.random.shuffle(years_data)
    return years_data

def H3N2_Npy_5folds(index,types):
    years_data = []
    for i in index:
        try:
            years_data.append(np.load("./training_data_5folds/IVA-CNN/{}_{}.npy".format(str(i),types)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data),axis=0)
    #np.random.shuffle(years_data)
    return years_data

def H3N2_CSV_Seq(years1,years2):
    temp = []
    for i in range(years1,years2+1):
        try:
            df = pd.read_csv('./sequence_nodup_single/data/{}.csv'.format(i),names=['seq', 'description'])
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans

def H3N2_CSV_Seq_5folds(index):
    temp = []
    for i in index:
        try:
            df = pd.read_csv('./sequence_5folds/data/{}.csv'.format(i),names=['seq', 'description'])
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans

def H3N2_CSV(path,years1,years2):
    temp = []
    for i in range(years1,years2+1):
        try:
            df = pd.read_csv(path+'{}.csv'.format(i))
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans

def H3N2_CSV_5folds(path,index):
    temp = []
    for i in index:
        try:
            df = pd.read_csv(path+'{}.csv'.format(i))
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans
    
def H3N2_CSV_Anti(years1,years2):
    temp = []
    for i in range(years1,years2+1):
        try:
            df = pd.read_csv('./antigen_nodup_single/data/{}.csv'.format(i))
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans

def H3N2_CSV_Anti_5folds(index):
    temp = []
    for i in index:
        try:
            df = pd.read_csv('./antigen_5folds/data/{}.csv'.format(i))
            temp.append(df)
        except:
            continue
    ans = pd.concat(temp,axis=0)
    return ans

def Shi_Lee_method(years1,years2):
    # Lee_method(1970,1980) 5 fold test
    # feature_type == 'Min-Shi Lee': 
    # don't need traing,  
    H3N2_Antigenic_dist = H3N2_CSV_Anti(years2,years2)
    H3N2_seq = H3N2_CSV_Seq(years2,years2)

    H3N2_num_mut_list = distance_mutation(H3N2_Antigenic_dist, H3N2_seq)    
    H3N2_Antigenic_dist_list = list(H3N2_Antigenic_dist['Distance'])
    conf_matrix = get_confusion_matrix(H3N2_num_mut_list, H3N2_Antigenic_dist_list, 'H3N2')        
    H3N2_acc = get_accuracy(conf_matrix)
    H3N2_pre = get_precision(conf_matrix)
    H3N2_rec = get_recall(conf_matrix)
    H3N2_f1 = get_f1score(conf_matrix)
    H3N2_mcc = get_mcc(conf_matrix)
    ans = '{},{},{},{},{},{}\n'.format(years1,H3N2_acc, H3N2_pre, H3N2_rec, H3N2_f1, H3N2_mcc) 
    ans= '{},'.format(years1) + '%.3f,%.3f,%.3f,%.3f,%.3f\n'% (H3N2_acc, H3N2_pre, H3N2_rec, H3N2_f1, H3N2_mcc)
    return ans

def Shi_Lee_5folds_method(index):
    # Lee_method(1970,1980) 5 fold test
    # feature_type == 'Min-Shi Lee': 
    # don't need traing,  
    H3N2_Antigenic_dist = H3N2_CSV_Anti_5folds(index)
    H3N2_seq = H3N2_CSV_Seq_5folds(index)

    H3N2_num_mut_list = distance_mutation(H3N2_Antigenic_dist, H3N2_seq)    
    H3N2_Antigenic_dist_list = list(H3N2_Antigenic_dist['Distance'])
    conf_matrix = get_confusion_matrix(H3N2_num_mut_list, H3N2_Antigenic_dist_list, 'H3N2')        
    H3N2_acc = get_accuracy(conf_matrix)
    H3N2_pre = get_precision(conf_matrix)
    H3N2_rec = get_recall(conf_matrix)
    H3N2_f1 = get_f1score(conf_matrix)
    H3N2_mcc = get_mcc(conf_matrix)
    ans = '{},{},{},{},{},{}\n'.format(str(index),H3N2_acc, H3N2_pre, H3N2_rec, H3N2_f1, H3N2_mcc) 
    ans= '{},'.format(str(index)) + '%.3f,%.3f,%.3f,%.3f,%.3f\n'% (H3N2_acc, H3N2_pre, H3N2_rec, H3N2_f1, H3N2_mcc)
    return ans

def Shi_Lee_test(years1,years2):
    # Lee_method(1970,1980) 5 fold test
    # feature_type == 'Min-Shi Lee': 
    # don't need traing,  
    # return pre_data, valid_true_data
    H3N2_Antigenic_dist = H3N2_CSV_Anti(years1,years2)
    H3N2_seq = H3N2_CSV_Seq(years1,years2)

    H3N2_num_mut_list = distance_mutation(H3N2_Antigenic_dist, H3N2_seq)    
    H3N2_Antigenic_dist_list = list(H3N2_Antigenic_dist['Distance'])
    conf_matrix = get_confusion_matrix_test(H3N2_num_mut_list, H3N2_Antigenic_dist_list, 'H3N2')        
    return conf_matrix
           
def Shi_Lee_test_ROC(years1,years2):
    # Lee_method(1970,1980) 5 fold test
    # feature_type == 'Min-Shi Lee': 
    # don't need traing,  
    # return pre_data, valid_true_data
    H3N2_Antigenic_dist = H3N2_CSV_Anti(years1,years2)
    H3N2_seq = H3N2_CSV_Seq(years1,years2)

    H3N2_num_mut_list = distance_mutation(H3N2_Antigenic_dist, H3N2_seq)    
    H3N2_Antigenic_dist_list = list(H3N2_Antigenic_dist['Distance'])
    conf_matrix = get_confusion_matrix_test(H3N2_num_mut_list, H3N2_Antigenic_dist_list, 'H3N2')        
    return conf_matrix


def Liao_method(years1,years2):
    # Yu-Chieh Liao :
    # Liao_method(1970,1972)
    path = './training_data/Yu-Chieh Liao/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any') # del the null value cow
                
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
                
    H3N2_data_train = H3N2_data_train.iloc[:,1:331]
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:331]

    H3N2_train_x = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_train_y = H3N2_data_train['label']
    H3N2_valid_x = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_valid_y = H3N2_data_valid['label']
    print('Yu-Chieh Liao method on H3N2 using svm:{}'.format(years2))
                #train_x, test_x, train_y, test_y = train_test_split_data(H3N2_feature, H3N2_label, 0.2)
    return svm_baseline(H3N2_train_x,H3N2_train_y,H3N2_valid_x,H3N2_valid_y)

def Liao_5folds_method(train,valid):
    # Yu-Chieh Liao :
    # Liao_method(1970,1972)
    path = './training_data_5folds/Yu-Chieh Liao/'
    H3N2_data_train = H3N2_CSV_5folds(path,train)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any') # del the null value cow
                
    H3N2_data_valid = H3N2_CSV_5folds(path,valid)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
                
    H3N2_data_train = H3N2_data_train.iloc[:,1:331]
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:331]

    H3N2_train_x = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_train_y = H3N2_data_train['label']
    H3N2_valid_x = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_valid_y = H3N2_data_valid['label']
    print('Yu-Chieh Liao method on H3N2 using svm:{}'.format(str(valid)))
                #train_x, test_x, train_y, test_y = train_test_split_data(H3N2_feature, H3N2_label, 0.2)
    return svm_baseline(H3N2_train_x,H3N2_train_y,H3N2_valid_x,H3N2_valid_y)

def Liao_test(years1,years2):
    # Yu-Chieh Liao :
    # Liao_method(1970,1972)
    path = './training_data/Yu-Chieh Liao/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any') # del the null value cow
                
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
                
    H3N2_data_train = H3N2_data_train.iloc[:,1:331]
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:331]

    H3N2_train_x = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_train_y = H3N2_data_train['label']
    H3N2_valid_x = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_valid_y = H3N2_data_valid['label']
    print('Yu-Chieh Liao method on H3N2 using svm:{}'.format(years2))
                #train_x, test_x, train_y, test_y = train_test_split_data(H3N2_feature, H3N2_label, 0.2)
    return svm_baseline_test(H3N2_train_x,H3N2_train_y,H3N2_valid_x,H3N2_valid_y)

def Liao_test_ROC(years1,years2):
    # Yu-Chieh Liao :
    # Liao_method(1970,1972)
    path = './training_data/Yu-Chieh Liao/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any') # del the null value cow
                
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
                
    H3N2_data_train = H3N2_data_train.iloc[:,1:331]
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:331]

    H3N2_train_x = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_train_y = H3N2_data_train['label']
    H3N2_valid_x = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_valid_y = H3N2_data_valid['label']
    print('Yu-Chieh Liao method on H3N2 using svm:{}'.format(years2))
                #train_x, test_x, train_y, test_y = train_test_split_data(H3N2_feature, H3N2_label, 0.2)
    return svm_baseline_test_ROC(H3N2_train_x,H3N2_train_y,H3N2_valid_x,H3N2_valid_y)


def William_Lee_method(years1,years2):
    # feature_type == 'William Lees':
    # Lee_method(1978,2013)
    path = './training_data/William Lees/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any')
    H3N2_data_train = H3N2_data_train.iloc[:,1:7]
    
    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_label_train = H3N2_data_train['label']
            
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:7]
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_label_valid = H3N2_data_valid['label']

    print('William Lees method on H3N2 using svm:{}'.format(years2))
    return svm_baseline(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def William_Lee_5folds_method(train,valid):
    # feature_type == 'William Lees':
    # Lee_method(1978,2013)
    path = './training_data_5folds/William Lees/'
    H3N2_data_train = H3N2_CSV_5folds(path,train)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any')
    H3N2_data_train = H3N2_data_train.iloc[:,1:7]
    
    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_label_train = H3N2_data_train['label']
            
    H3N2_data_valid = H3N2_CSV_5folds(path,valid)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:7]
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_label_valid = H3N2_data_valid['label']

    print('William Lees method on H3N2 using svm:{}'.format(str(valid)))
    return svm_baseline(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def William_Lee_test(years1,years2):
    # feature_type == 'William Lees':
    # Lee_method(1978,2013)
    path = './training_data/William Lees/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any')
    H3N2_data_train = H3N2_data_train.iloc[:,1:7]
    
    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_label_train = H3N2_data_train['label']
            
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:7]
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_label_valid = H3N2_data_valid['label']

    print('William Lees method on H3N2 using svm:{}'.format(years2))
    return svm_baseline_test(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def William_Lee_test_ROC(years1,years2):
    # feature_type == 'William Lees':
    # Lee_method(1978,2013)
    path = './training_data/William Lees/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any')
    H3N2_data_train = H3N2_data_train.iloc[:,1:7]
    
    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label']
    H3N2_label_train = H3N2_data_train['label']
            
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:7]
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label']
    H3N2_label_valid = H3N2_data_valid['label']

    print('William Lees method on H3N2 using svm:{}'.format(years2))
    return svm_baseline_test_ROC(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )



def Peng_method(years1,years2):
    #'Peng Yousong': Peng Yousong , 5 fold test
    # Peng_method(1969,1973)
    path = './training_data/Peng Yousong/'
    H3N2_data = H3N2_CSV(path,years1,years2)
    H3N2_data = H3N2_data.dropna(axis=0, how='any')
    H3N2_data = H3N2_data.iloc[:,1:12]
    
    H3N2_feature = H3N2_data.iloc[:, H3N2_data.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label = H3N2_data['label'].reset_index().iloc[:,1:]
    
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:12]
    
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label_valid = H3N2_data_valid['label'].reset_index().iloc[:,1:]
    H3N2_label_valid = np.array(list(H3N2_data_valid['label']))

    print('Peng Yousong method on H3N2 using naive bayes:{}'.format(years2+1))
    return bayes_cross_validation(H3N2_feature, H3N2_label, H3N2_feature_valid, H3N2_label_valid)

def Peng_5folds_method(train,valid):
    #'Peng Yousong': Peng Yousong , 5 fold test
    # Peng_method(1969,1973)
    path = './training_data_5folds/Peng Yousong/'
    H3N2_data = H3N2_CSV_5folds(path,train)
    H3N2_data = H3N2_data.dropna(axis=0, how='any')
    H3N2_data = H3N2_data.iloc[:,1:12]
    
    H3N2_feature = H3N2_data.iloc[:, H3N2_data.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label = H3N2_data['label'].reset_index().iloc[:,1:]
    
    H3N2_data_valid = H3N2_CSV_5folds(path,valid)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:12]
    
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label_valid = H3N2_data_valid['label'].reset_index().iloc[:,1:]
    H3N2_label_valid = np.array(list(H3N2_data_valid['label']))

    print('Peng Yousong method on H3N2 using naive bayes:{}'.format(str(valid)))
    return bayes_cross_validation(H3N2_feature, H3N2_label, H3N2_feature_valid, H3N2_label_valid)

def Peng_test(years1,years2):
    #'Peng Yousong': Peng Yousong , 5 fold test
    # Peng_method(1969,1973)
    path = './training_data/Peng Yousong/'
    H3N2_data = H3N2_CSV(path,years1,years2)
    H3N2_data = H3N2_data.dropna(axis=0, how='any')
    H3N2_data = H3N2_data.iloc[:,1:12]
    
    H3N2_feature = H3N2_data.iloc[:, H3N2_data.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label = H3N2_data['label'].reset_index().iloc[:,1:]
    
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:12]
    
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label_valid = H3N2_data_valid['label'].reset_index().iloc[:,1:]
    H3N2_label_valid = np.array(list(H3N2_data_valid['label']))

    return bayes_cross_validation_test(H3N2_feature, H3N2_label, H3N2_feature_valid, H3N2_label_valid)
    
def Peng_test_ROC(years1,years2):
    #'Peng Yousong': Peng Yousong , 5 fold test
    # Peng_method(1969,1973)
    path = './training_data/Peng Yousong/'
    H3N2_data = H3N2_CSV(path,years1,years2)
    H3N2_data = H3N2_data.dropna(axis=0, how='any')
    H3N2_data = H3N2_data.iloc[:,1:12]
    
    H3N2_feature = H3N2_data.iloc[:, H3N2_data.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label = H3N2_data['label'].reset_index().iloc[:,1:]
    
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any')
    H3N2_data_valid = H3N2_data_valid.iloc[:,1:12]
    
    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index().iloc[:,1:]
    H3N2_label_valid = H3N2_data_valid['label'].reset_index().iloc[:,1:]
    H3N2_label_valid = np.array(list(H3N2_data_valid['label']))

    return bayes_cross_validation_test_ROC(H3N2_feature, H3N2_label, H3N2_feature_valid, H3N2_label_valid)
    


def Yao_method(years1,years2):
    # 'Yuhua Yao':
    # Yao_method(1979,1983)   
    path = './training_data/Yuhua Yao/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any').reset_index()
    
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any').reset_index()
                
    H3N2_data_train = H3N2_data_train.iloc[:,2:332]
    H3N2_data_valid = H3N2_data_valid.iloc[:,2:332]

    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label'].reset_index()
    H3N2_label_train = H3N2_data_train['label'] #.reset_index()

    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index()
    H3N2_label_valid = H3N2_data_valid['label'] #.reset_index()
    print('Yuhua Yao method on H3N2 using random forest:{}'.format(years2+1))
    #randomforest_cross_validation(H3N2_feature, H3N2_label)
    return rf_baseline(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def Yao_5folds_method(train,valid):
    # 'Yuhua Yao':
    # Yao_method(1979,1983)   
    path = './training_data_5folds/Yuhua Yao/'
    H3N2_data_train = H3N2_CSV_5folds(path,train)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any').reset_index()
    
    H3N2_data_valid = H3N2_CSV_5folds(path,valid)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any').reset_index()
                
    H3N2_data_train = H3N2_data_train.iloc[:,2:332]
    H3N2_data_valid = H3N2_data_valid.iloc[:,2:332]

    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label'].reset_index()
    H3N2_label_train = H3N2_data_train['label'] #.reset_index()

    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index()
    H3N2_label_valid = H3N2_data_valid['label'] #.reset_index()
    #randomforest_cross_validation(H3N2_feature, H3N2_label)
    return rf_baseline(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def Yao_test(years1,years2):
    # 'Yuhua Yao':
    # Yao_method(1979,1983)   
    path = './training_data/Yuhua Yao/'
    H3N2_data_train = H3N2_CSV(path,years1,years2)
    H3N2_data_train = H3N2_data_train.dropna(axis=0, how='any').reset_index()
    
    H3N2_data_valid = H3N2_CSV(path,years2+1,years2+1)
    H3N2_data_valid = H3N2_data_valid.dropna(axis=0, how='any').reset_index()
                
    H3N2_data_train = H3N2_data_train.iloc[:,2:332]
    H3N2_data_valid = H3N2_data_valid.iloc[:,2:332]

    H3N2_feature_train = H3N2_data_train.iloc[:, H3N2_data_train.columns!='label'].reset_index()
    H3N2_label_train = H3N2_data_train['label'] #.reset_index()

    H3N2_feature_valid = H3N2_data_valid.iloc[:, H3N2_data_valid.columns!='label'].reset_index()
    H3N2_label_valid = H3N2_data_valid['label'] #.reset_index()
    print('Yuhua Yao method on H3N2 using random forest:{}'.format(years2+1))
    #randomforest_cross_validation(H3N2_feature, H3N2_label)
    return rf_baseline_test(H3N2_feature_train, H3N2_label_train, H3N2_feature_valid , H3N2_label_valid )

def IAV_CNN_method(years1,years2):
    
    parameters = { 
      # 'rf', lr', 'knn', 'svm', 'cnn'
      #'model': model,
      # Number of hidden units in the encoder
      'hidden_size': 128,
      # Droprate (applied at input)
      'dropout_p': 0.5,
      # Note, no learning rate decay implemented
      'learning_rate': 0.001,
      # Size of mini batch
      'batch_size': 32,
      # Number of training iterations
      'num_of_epochs': 100
    }
    train_x = H3N2_Npy(years1, years2 ,'x')
    train_y = H3N2_Npy(years1, years2 ,'y')

    test_x = H3N2_Npy(years2+1, years2+1 ,'x')
    test_y = H3N2_Npy(years2+1, years2+1 ,'y')

    train_x = np.reshape(train_x, (np.array(train_x).shape[0], 1, np.array(train_x).shape[1], np.array(train_x).shape[2]))
    test_x =  np.reshape(test_x, (np.array(test_x).shape[0], 1, np.array(test_x).shape[1], np.array(test_x).shape[2]))
    #print(np.array(train_x).shape)
    #print(np.array(test_x).shape)
            
    train_x = torch.tensor(train_x, dtype=torch.float32).cuda()
    train_y = torch.tensor(train_y, dtype=torch.int64).cuda()
    test_x = torch.tensor(test_x, dtype=torch.float32).cuda()
    test_y = torch.tensor(test_y, dtype=torch.int64).cuda()

    #net = CNN_H3N2(1, 128, 2, 2)
    net = CNN_H3N2()
    net.cuda()
    return train_cnn(net, parameters['num_of_epochs'], parameters['learning_rate'], parameters['batch_size'], train_x, train_y, test_x, test_y)     

def IAV_CNN_test(years1,years2):
    
    parameters = { 
      # 'rf', lr', 'knn', 'svm', 'cnn'
      #'model': model,
      # Number of hidden units in the encoder
      'hidden_size': 128,
      # Droprate (applied at input)
      'dropout_p': 0.5,
      # Note, no learning rate decay implemented
      'learning_rate': 0.001,
      # Size of mini batch
      'batch_size': 32,
      # Number of training iterations
      'num_of_epochs': 100
    }
    train_x = H3N2_Npy(years1, years2 ,'x')
    train_y = H3N2_Npy(years1, years2 ,'y')

    test_x = H3N2_Npy(years2+1, years2+1 ,'x')
    test_y = H3N2_Npy(years2+1, years2+1 ,'y')

    train_x = np.reshape(train_x, (np.array(train_x).shape[0], 1, np.array(train_x).shape[1], np.array(train_x).shape[2]))
    test_x =  np.reshape(test_x, (np.array(test_x).shape[0], 1, np.array(test_x).shape[1], np.array(test_x).shape[2]))
    #print(np.array(train_x).shape)
    #print(np.array(test_x).shape)
            
    train_x = torch.tensor(train_x, dtype=torch.float32).cuda()
    train_y = torch.tensor(train_y, dtype=torch.int64).cuda()
    test_x = torch.tensor(test_x, dtype=torch.float32).cuda()
    test_y = torch.tensor(test_y, dtype=torch.int64).cuda()

    #net = CNN_H3N2(1, 128, 2, 2)
    net = CNN_H3N2()
    net.cuda()
    return test_cnn(net, parameters['num_of_epochs'], parameters['learning_rate'], parameters['batch_size'], train_x, train_y, test_x, test_y)     

def IAV_CNN_test_ROC(years1,years2):
    
    parameters = { 
      # 'rf', lr', 'knn', 'svm', 'cnn'
      #'model': model,
      # Number of hidden units in the encoder
      'hidden_size': 128,
      # Droprate (applied at input)
      'dropout_p': 0.5,
      # Note, no learning rate decay implemented
      'learning_rate': 0.001,
      # Size of mini batch
      'batch_size': 32,
      # Number of training iterations
      'num_of_epochs': 100
    }
    train_x = H3N2_Npy(years1, years2 ,'x')
    train_y = H3N2_Npy(years1, years2 ,'y')

    test_x = H3N2_Npy(years2+1, years2+1 ,'x')
    test_y = H3N2_Npy(years2+1, years2+1 ,'y')

    train_x = np.reshape(train_x, (np.array(train_x).shape[0], 1, np.array(train_x).shape[1], np.array(train_x).shape[2]))
    test_x =  np.reshape(test_x, (np.array(test_x).shape[0], 1, np.array(test_x).shape[1], np.array(test_x).shape[2]))
    #print(np.array(train_x).shape)
    #print(np.array(test_x).shape)
            
    train_x = torch.tensor(train_x, dtype=torch.float32).cuda()
    train_y = torch.tensor(train_y, dtype=torch.int64).cuda()
    test_x = torch.tensor(test_x, dtype=torch.float32).cuda()
    test_y = torch.tensor(test_y, dtype=torch.int64).cuda()

    #net = CNN_H3N2(1, 128, 2, 2)
    net = CNN_H3N2()
    net.cuda()
    return test_cnn_ROC(net, parameters['num_of_epochs'], parameters['learning_rate'], parameters['batch_size'], train_x, train_y, test_x, test_y)     



def IAV_CNN_5folds_method(train,valid):
    
    parameters = { 
      # 'rf', lr', 'knn', 'svm', 'cnn'
      #'model': model,
      # Number of hidden units in the encoder
      'hidden_size': 128,
      # Droprate (applied at input)
      'dropout_p': 0.5,
      # Note, no learning rate decay implemented
      'learning_rate': 0.001,
      # Size of mini batch
      'batch_size': 32,
      # Number of training iterations
      'num_of_epochs': 100
    }
    train_x = H3N2_Npy_5folds(train ,'x')
    train_y = H3N2_Npy_5folds(train ,'y')

    test_x = H3N2_Npy_5folds(valid ,'x')
    test_y = H3N2_Npy_5folds(valid ,'y')

    train_x = np.reshape(train_x, (np.array(train_x).shape[0], 1, np.array(train_x).shape[1], np.array(train_x).shape[2]))
    test_x =  np.reshape(test_x, (np.array(test_x).shape[0], 1, np.array(test_x).shape[1], np.array(test_x).shape[2]))
    #print(np.array(train_x).shape)
    #print(np.array(test_x).shape)
            
    train_x = torch.tensor(train_x, dtype=torch.float32).cuda()
    train_y = torch.tensor(train_y, dtype=torch.int64).cuda()
    test_x = torch.tensor(test_x, dtype=torch.float32).cuda()
    test_y = torch.tensor(test_y, dtype=torch.int64).cuda()

    #net = CNN_H3N2(1, 128, 2, 2)
    net = CNN_H3N2()
    net.cuda()
    return train_cnn(net, parameters['num_of_epochs'], parameters['learning_rate'], parameters['batch_size'], train_x, train_y, test_x, test_y)     

    