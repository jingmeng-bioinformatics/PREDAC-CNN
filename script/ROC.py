from sklearn.metrics import roc_curve, auc, roc_auc_score
import numpy as np
import matplotlib.pyplot as plt
from utils import *

def draw_ROC(years1, years2):
    methods = ['H3N2-CNN',  'H3N2-CNN_only_single', 'Yu-Chieh Liao', 'Du_Xiangjun', 'William Lees', 'IAV-CNN']  # 'Min-Shi Lee',
    methods_name = ['CNNantigenic',  'CNNantigenic-single', 'Liao', 'PREDAC', 'Lees', 'IAV-CNN']  # 'Min-Shi Lee',
    colors = ['lightskyblue', 'darkorange', 'darkgrey', 'gold', 'cornflowerblue', 'darkblue', 'yellowgreen']
    plt.figure(figsize=(7, 6))
    for j in range(len(methods)):
        years_data = []
        for i in range(years1, years2 + 1):
            years_data.append(np.load("./ROC/{}/{}.npy".format(methods[j], i), allow_pickle=True))
        a = np.concatenate(tuple(years_data), axis=1)
        if(methods[j]=='Min-Shi Lee' or methods[j]=='Du_Xiangjun'):
            auc = roc_auc_score(a[1, :], a[0, :])
            fpr, tpr, thresholds = roc_curve(a[1, :], a[0, :])
        else:
            auc = roc_auc_score(a[1, :, 0], a[0, :, 0])
            fpr, tpr, thresholds = roc_curve(a[1, :, 0], a[0, :, 0])

        plt.plot(fpr, tpr, color=colors[j], label='{} ({:.2f})'.format(methods_name[j], auc))
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('{}-{}'.format(years1,years2))
        plt.legend(loc="lower right")
    plt.savefig('./ROC/years/{}_{}.png'.format(years1, years2), dpi=500)
    plt.clf()

def Score(years1, years2):
    methods = ['H3N2-CNN',  'H3N2-CNN_only_single', 'Yu-Chieh Liao', 'Du_Xiangjun', 'William Lees', 'IAV-CNN']  # 'Min-Shi Lee',
    methods_name = ['CNNantigenic',  'CNNantigenic-single', 'Liao', 'PREDAC', 'Lees', 'IAV-CNN']  # 'Min-Shi Lee',
    f = open("./Judge_score/years_all_model_{}-{}.csv".format(years1, years2), "w")
    f.write("Methods,precision, recall, fscore, mcc, val_acc\n")
    for j in range(len(methods)):
        years_data = []
        for i in range(years1, years2 + 1):
            years_data.append(np.load("./ROC/{}/{}.npy".format(methods[j], i), allow_pickle=True))
        a = np.concatenate(tuple(years_data), axis=1)
        if(methods[j] != 'Min-Shi Lee' and methods[j] != 'Du_Xiangjun'):
            temp2 = Format_Y2(a[0], a[1])
        else:
            temp2 = a
        precision, recall, fscore, mcc, val_acc = evaluate(temp2[0], temp2[1])
        f.write("{},{}, {}, {}, {}, {}\n".format(methods_name[j],precision, recall, fscore, mcc, val_acc))
    f.close()

def TPFP():
    methods = ['Du_Xiangjun']
    for j in range(len(methods)):
        for i in range(2006, 2021):
            a=np.load("./ROC/{}/{}.npy".format(methods[j], i))
            if(methods[j] != 'Min-Shi Lee' and methods[j] != 'Du_Xiangjun'):
                temp2 = Format_Y2(a[0], a[1])
            else:
                temp2 = a
            precision, recall, fscore, mcc, val_acc = evaluate(temp2[0], temp2[1])
            print(val_acc)
        print("~~~~~~~")
TPFP()
'''
for i in range(2006,2007):
    draw_ROC(i, i+2)
    Score(i, i+2)

for i in range(2009, 2019, 4):
    print(i)
    draw_ROC(i, i+3)
    Score(i, i+3)
'''

