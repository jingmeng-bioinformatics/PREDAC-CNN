import os
#years train
choice = '5folds' #'years'
if(choice =='years'):
    for i in range(2000,2020):
        os.system("perl modelGenerator.pl epitope_H3N2 aaIndex receptor_H3N2 ..\\training_data\Du_Xiangjun\H3N2_seq_1968_{} ..\\training_data\Du_Xiangjun\seqCompare_H3N2_1968_{} .\\save_model\moel_1968_{}".format(i,i,i))

if(choice =='5folds'):
    for i in range(5):
        os.system("perl modelGenerator.pl epitope_H3N2 aaIndex receptor_H3N2 ..\\training_data_5folds\Du_Xiangjun\H3N2_seq_train_{} ..\\training_data_5folds\Du_Xiangjun\seqCompare_H3N2_train_{} .\\save_model\moel_train_{}".format(i,i,i))

