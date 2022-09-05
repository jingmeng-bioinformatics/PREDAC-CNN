import os, sys
import pandas as pd
import numpy as np
import random
import torch
import warnings
import math


from Bio import SeqIO
from model import fasta_to_csv
from model import calculate_label
from model import generate_feature
from model import strain_selection
from model import replace_uncertain_amino_acids
from model import train_test_split_data

warnings.filterwarnings('ignore')

#William method
def William(Antigenic_dist,seq,save_csv):
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq,names=['seq', 'description'])
    H3N2_new_epitope_a = [71, 72, 98, 122, 124, 126, 127, 130, 131, 132, 133, 135, 137, 138, 140, 141, 142, 143, 144, 145, 146, 148, 149,
                          150, 151, 152, 168, 255]
    H3N2_new_epitope_b = [128, 129, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 
                          197, 198, 199]
    H3N2_new_epitope_c = [41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 271, 272, 273, 274, 275, 276, 278, 279, 280, 282, 
                          284, 285, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 307, 
                          308, 309, 310, 311, 312, 313, 314]
    H3N2_new_epitope_d = [95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 107, 117, 118, 120, 121, 166, 167, 169, 170, 171, 172, 173, 174, 
                          175, 176, 177, 178, 179, 180, 182, 183, 184, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 
                          214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 
                          236, 238, 239, 240, 242, 243, 244, 245, 246, 247, 248, 257, 258]
    H3N2_new_epitope_e = [56, 57, 58, 59, 60, 62, 63, 64, 65, 67, 68, 69, 70, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87,
                          88, 89, 90, 91, 92, 93, 94, 109, 110, 111, 112, 113, 114, 115, 119, 259, 260, 261, 262, 263, 264, 265, 267, 268,
                          268, 270]
    H3N2_new_epitopes = {'new_epitope_a': H3N2_new_epitope_a, 'new_epitope_b': H3N2_new_epitope_b, 'new_epitope_c': H3N2_new_epitope_c, 
                    'new_epitope_d': H3N2_new_epitope_d, 'new_epitope_e': H3N2_new_epitope_e}
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq, names=['seq', 'description'])
    print("generate_feature({},{},{})...".format(Antigenic_dist,seq,save_csv))
    H3N2_epitope_data = generate_feature(H3N2_new_epitopes, H3N2_Antigenic_dist, H3N2_seq)
    H3N2_epitope_data.to_csv(save_csv)
    #'training/William Lees/H3N2_new_epitope_data.csv'


#
#
##0 is H1N1, 1 is H3N2, 2 is H5N1, map antigenic site across different subtypes
##embedding = pd.read_csv('sequence/embedding.csv', header=None)
def map_antigenic_site(embedding_table, input_site, subtype):
    input_index = []
    output_index = []
    if subtype == 1:
        input_index = list(embedding_table[8])
    for i in range(0, len(input_site)):
        print("---{:.2f}---".format(100*i/len(input_site)),end="\r")
        for j in range(embedding_table.shape[0]):
            result = []
            if input_index[j] == str(input_site[i]):
                result = int(embedding_table[12][j])
                break                
        output_index.append(result)
    return output_index


#Peng Yousong's work
def Peng(Antigenic_dist,seq,save_csv):
    #divide residue sites based on ten regional bands and generate training data with new features
    #H3N2
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq,names=['seq', 'description'])
    H3N2_regional_1 = [129,130,131,155,156,157,158,159,160,196]
    H3N2_regional_2 = [126,127,128,132,133,134,135,136,153,162,163,165,188,189,190,192,193,194,197,198,199,201,248]
    H3N2_regional_3 = [74,100,101,122,123,124,125,137,138,140,141,142,143,144,145,147,149,150,166,167,168,186,187,
                       211,212,213,214,215,216,217,219,222,223,224,225,226,227,233,234,244,255,257]
    H3N2_regional_4 = [63,65,75,77,78,79,80,81,82,93,94,95,96,102,103,104,105,106,109,119,120,121,169,171,172,174,
                       175,207,208,209,210,236,238,239,240,242,259,260]
    H3N2_regional_5 = [57,58,59,60,62,83,85,89,90,91,92,110,114,173,261,262,263,264,267,269,271]
    H3N2_regional_6 = [49,50,53,54,55,56,272,273,274,275,276,277,278,279,280,284,285,298,299,300,301]
    H3N2_regional_7 = [44,45,46,47,48,289,290,291,292,293,296,297,307,308,310,311,312]
    H3N2_regional_8 = [40,41,313,315]
    H3N2_regional_9 = [22,23,24,25,27,29,31,32,33,34,35,37,38,39,318]
    H3N2_regional_10 = [1,2,3,4,5,6,7,8,9,10,12,14,18,20,21,321,323,324,325,326,327,328]
    H3N2_regional_band = {'regional_1':H3N2_regional_1, 'regional_2':H3N2_regional_2, 'regional_3':H3N2_regional_3, 'regional_4':H3N2_regional_4,
                          'regional_5':H3N2_regional_5, 'regional_6':H3N2_regional_6, 'regional_7':H3N2_regional_7, 'regional_8':H3N2_regional_8, 
                          'regional_9':H3N2_regional_9, 'regional_10':H3N2_regional_10}

    print("Peng generate_feature({},{},{})...".format(Antigenic_dist,seq,save_csv))
    H3N2_regional_data = generate_feature(H3N2_regional_band, H3N2_Antigenic_dist, H3N2_seq)
    # 'training/Peng Yousong/H3N2_regional_data.csv'
    H3N2_regional_data.to_csv(save_csv)

def Liao_feature_engineering(distance_input, seq_input, subtype):
    group1 = ['A', 'I', 'L', 'M', 'P', 'V']
    group2 = ['F', 'W', 'Y']       
    group3 = ['N', 'Q', 'S', 'T']
    group4 = ['D', 'E', 'H', 'K', 'R']
    group5 = ['C']
    group6 = ['G']
    distance_label = calculate_label(distance_input)
    label = {'label': distance_label}
    label = pd.DataFrame(label)
    
    index = pd.Series(np.arange(distance_input.shape[0]))
    length = len(seq_input['seq'].iloc[0])

    if subtype == 0:
        columns = list(range(1, 328, 1))
    elif subtype == 1:
        columns = list(range(1, 330, 1))
    elif subtype == 2:
        columns = list(range(1, 321, 1))
    for col in range(len(columns)):
        columns[col] = str(columns[col])
    Mut_feature = pd.DataFrame(index=index, columns=columns)
    
    for i in range(0, distance_input.shape[0]):
        print("---{}---".format(100*i/distance_input.shape[0]))
        strain_1 = []
        strain_2 = []
        for j in range(0, seq_input.shape[0]):
            if seq_input['description'].iloc[j] == distance_input['Strain1'].iloc[i]:
                strain_1 = seq_input['seq'].iloc[j].upper()
            if seq_input['description'].iloc[j] == distance_input['Strain2'].iloc[i]:
                strain_2 = seq_input['seq'].iloc[j].upper()
                
        for a in range(0, length):
            if strain_1[a] in group1 and strain_2[a] in group1:
                Mut_feature.iloc[i][a] = 0
            elif strain_1[a] in group2 and strain_2[a] in group2:
                Mut_feature.iloc[i][a] = 0
            elif strain_1[a] in group3 and strain_2[a] in group3:
                Mut_feature.iloc[i][a] = 0
            elif strain_1[a] in group4 and strain_2[a] in group4:
                Mut_feature.iloc[i][a] = 0
            elif strain_1[a] in group5 and strain_2[a] in group5:
                Mut_feature.iloc[i][a] = 0
            elif strain_1[a] in group6 and strain_2[a] in group6:
                Mut_feature.iloc[i][a] = 0
            else:
               Mut_feature.iloc[i][a] = 1
    Mut_feature = Mut_feature.join(label)   
    return(Mut_feature)


#Yu-Chieh Liao method (Bioinformatics models for predicting antigenic variants of influenza A/H3N2 virus)
def Liao(Antigenic_dist,seq,save_csv):
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq, names=['seq', 'description'])
    print("Liao_feature_H3N2 ({},{},{})...".format(Antigenic_dist,seq,save_csv))
    Liao_feature_H3N2 = Liao_feature_engineering(H3N2_Antigenic_dist, H3N2_seq, 1)
    Liao_feature_H3N2.to_csv(save_csv)


def Yao_feature_engineering(distance_input, seq_input, subtype):
    Yao_embedding = pd.read_csv('./input_map/NIEK910102_matrix.csv')
    Yao_pair_aa = []
    Yao_pair_aa_value = []
    for i in range(1, Yao_embedding.shape[0]): #start from index 1
        print("---{:.2f}---".format(100*i/Yao_embedding.shape[0]),end="\r")
        aa_pair = 0
        aa_value = 0
        for j in range(1, Yao_embedding.shape[1]):
            aa_pair = Yao_embedding.iloc[i][0] + Yao_embedding.iloc[0][j]
            aa_value = Yao_embedding.iloc[i][j]
            if not math.isnan(float(aa_value)):
                Yao_pair_aa.append(aa_pair)
                Yao_pair_aa_value.append(aa_value)
    Yao_embedding_table = pd.DataFrame({'aa_pair': Yao_pair_aa, 'aa_value': Yao_pair_aa_value})
    distance_label = calculate_label(distance_input)
    label = {'label': distance_label}
    label = pd.DataFrame(label)
    
    index = pd.Series(np.arange(distance_input.shape[0]))
    length = len(seq_input['seq'].iloc[0])    
    if subtype == 0:
        columns = list(range(1, 328, 1))
    elif subtype == 1:
        columns = list(range(1, 330, 1))
    elif subtype == 2:
        columns = list(range(1, 321, 1))
    for col in range(len(columns)):
        columns[col] = str(columns[col])
    Mut_feature = pd.DataFrame(index=index, columns=columns)

    for i in range(0, distance_input.shape[0]):
        print("---{:.2f}---".format(100*i/distance_input.shape[0]),end="\r")
        strain_1 = []
        strain_2 = []
        for j in range(0, seq_input.shape[0]):
            if seq_input['description'].iloc[j] == distance_input['Strain1'].iloc[i]:
                strain_1 = seq_input['seq'].iloc[j].upper()
            if seq_input['description'].iloc[j] == distance_input['Strain2'].iloc[i]:
                strain_2 = seq_input['seq'].iloc[j].upper()
        for a in range(0, length):
            aa_pair_1 = strain_1[a] + strain_2[a]
            aa_pair_2 = strain_2[a] + strain_1[a]
            for b in range(0, Yao_embedding_table.shape[0]):
                if strain_1[a] == '-' or strain_2[a] == '-':
                    Mut_feature.iloc[i][a] = 0
                    break
                elif aa_pair_1 == Yao_embedding_table['aa_pair'].iloc[b] or aa_pair_2 == Yao_embedding_table['aa_pair'].iloc[b]:
                    Mut_feature.iloc[i][a] = Yao_embedding_table['aa_value'].iloc[b]
                    break

    Mut_feature = Mut_feature.join(label)
    return(Mut_feature)

##Yuhua Yao method (Predicting influenza antigenicity from Hemagglutintin sequence data based on a joint random forest method)
def Yao(Antigenic_dist,seq,save_csv):
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq, names=['seq', 'description'])
    print("Yao_feature_engineering({},{},{})...".format(Antigenic_dist,seq,save_csv))
    H3N2_seq = replace_uncertain_amino_acids(H3N2_seq)
    Yao_feature_H3N2 = Yao_feature_engineering(H3N2_Antigenic_dist, H3N2_seq, 1)
    Yao_feature_H3N2.to_csv(save_csv)



###################################################################################################################################
def cnn_training_data(Antigenic_dist, seq):
    raw_data = strain_selection(Antigenic_dist, seq)
    #replace unambiguous with substitutions
    Btworandom = 'DN'
    Jtworandom = 'IL'
    Ztworandom = 'EQ'
    Xallrandom = 'ACDEFGHIKLMNPQRSTVWY'
    for i in range(0, 2):
        for j in range(0, len(raw_data[0])):
            print("---{:.2f}---".format(j/len(raw_data[0])),end="\r")
            seq = raw_data[i][j]
            seq = seq.replace('B', random.choice(Btworandom))
            seq = seq.replace('J', random.choice(Jtworandom))
            seq = seq.replace('Z', random.choice(Ztworandom))
            seq = seq.replace('X', random.choice(Xallrandom))
            raw_data[i][j] = seq
            
    #embedding with ProVect    
    df = pd.read_csv('./input_map/protVec_100d_3grams.csv', delimiter = '\t')
    trigrams = list(df['words'])
    trigram_to_idx = {trigram: i for i, trigram in enumerate(trigrams)}
    trigram_vecs = df.loc[:, df.columns != 'words'].values

    feature = []
    label = raw_data[2]
    for i in range(0, len(raw_data[0])):
        print("---{:.2f}---".format(i/len(raw_data[0])),end="\r")
        trigram1 = []
        trigram2 = []
        strain_embedding = []
        seq1 = raw_data[0][i]
        seq2 = raw_data[1][i]
    
        for j in range(0, len(raw_data[0][0])-2):
            trigram1 = seq1[j:j+3]
            if trigram1[0] == '-' or trigram1[1] == '-' or trigram1[2] == '-':
                tri1_embedding = trigram_vecs[trigram_to_idx['<unk>']]
            else:
                tri1_embedding = trigram_vecs[trigram_to_idx[trigram1]]
        
            trigram2 = seq2[j:j+3]
            if trigram2[0] == '-' or trigram2[1] == '-' or trigram2[2] == '-':
                tri2_embedding = trigram_vecs[trigram_to_idx['<unk>']]
            else:
                tri2_embedding = trigram_vecs[trigram_to_idx[trigram2]]
        
            tri_embedding = tri1_embedding - tri2_embedding
            strain_embedding.append(tri_embedding)
        feature.append(strain_embedding)
    return feature, label

def IAV(Antigenic_dist,seq,save_x,save_y):
    H3N2_Antigenic_dist = pd.read_csv(Antigenic_dist)
    H3N2_seq = pd.read_csv(seq,names=['seq', 'description'])
    print("IAV_cnn_training_data({},{},{})...".format(Antigenic_dist,seq,save_x))
    ff,ll = cnn_training_data(H3N2_Antigenic_dist,H3N2_seq)
    ff = np.array(ff)
    ll = np.array(ll)
    np.save(save_x,ff)
    np.save(save_y,ll)

def Du(Relation,seq,save,year):
    H3N2_Relation = pd.read_csv(Relation)
    H3N2_seq = pd.read_csv(seq,names=['seq', 'description'])
    f = open(save+"H3N2_seq_{}".format(year),"w")
    for i in range(len(H3N2_seq)):
        f.write(">"+H3N2_seq['description'][i]+"\n")
        f.write(H3N2_seq['seq'][i]+"\n")
    f.close()
    
    f = open(save+"seqCompare_H3N2_{}".format(year),"w")
    for i in range(len(H3N2_Relation)):
        f.write(H3N2_Relation["Pair1"][i])
        f.write("\t"+H3N2_Relation["Pair2"][i])
        f.write("\t"+str(H3N2_Relation["Class"][i]))
        seq1 = H3N2_Relation["Seq1"][i]
        seq2 = H3N2_Relation["Seq2"][i]
        for j in range(329):
            if(seq1[j]=='-' or seq2[j]=='-'):
                continue
            elif(seq1[j]==seq2[j]):
                continue
            else:
                f.write("\t"+seq1[j]+str(j+1)+seq2[j])
        f.write("\n")
    f.close()
    print("Du Success")
    
def Du_train():
    #generate the Du train data
    save = './training_data/Du_Xiangjun/'
    for years in range(2000,2020):
        f1 = open(save+"H3N2_seq_{}_{}".format(1968,years),"w")
        f2 = open(save+"seqCompare_H3N2_{}_{}".format(1968,years),"w")
        for year in range(1968,years+1):
            try:
                Relation = "./CNN_nodup_years/{}.csv".format(year)
                H3N2_Relation = pd.read_csv(Relation)
                seq = './sequence_nodup_single/data/{}.csv'.format(year)
                H3N2_seq = pd.read_csv(seq,names=['seq', 'description'])

                for i in range(len(H3N2_seq)):
                    f1.write(">"+H3N2_seq['description'][i]+"\n")
                    f1.write(H3N2_seq['seq'][i]+"\n")

                for i in range(len(H3N2_Relation)):
                    f2.write(H3N2_Relation["Pair1"][i])
                    f2.write("\t"+H3N2_Relation["Pair2"][i])
                    f2.write("\t"+str(H3N2_Relation["Class"][i]))
                    seq1 = H3N2_Relation["Seq1"][i]
                    seq2 = H3N2_Relation["Seq2"][i]
                    for j in range(329):
                        if(seq1[j]=='-' or seq2[j]=='-'):
                            continue
                        elif(seq1[j]==seq2[j]):
                            continue
                        else:
                            f2.write("\t"+seq1[j]+str(j+1)+seq2[j])
                    f2.write("\n")
            except:
                continue
        f1.close()
        f2.close()            
    print("Du Success")


def main():
    for i in range(2000,2021):
        try:
            print("Year: {} ".format(i))
            Antigenic_dist = './antigen_nodup_single/data/{}.csv'.format(i)
            seq = './sequence_nodup_single/data/{}.csv'.format(i)
            '''
            William(Antigenic_dist,seq,'./training_data/William Lees/{}.csv'.format(str(i)))
            Peng(Antigenic_dist,seq,'./training_data/Peng Yousong/{}.csv'.format(str(i)))
            Liao(Antigenic_dist,seq,'./training_data/Yu-Chieh Liao/{}.csv'.format(str(i)))
            Yao(Antigenic_dist,seq,'./training_data/Yuhua Yao/{}.csv'.format(str(i)))
            IAV(Antigenic_dist,seq,'./training_data/IVA-CNN/{}_x.npy'.format(str(i)),'./training_data/IVA-CNN/{}_y.npy'.format(str(i)))
            '''
            Relation = "./CNN_nodup_years/{}.csv".format(i)
            Du(Relation,seq,'./training_data/Du_Xiangjun/',i)
        except:
            continue
Du_train()




















