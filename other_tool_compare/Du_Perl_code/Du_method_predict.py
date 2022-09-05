import os
choice = '5folds' #'years'
if(choice =='years'):
    for i in range(2000,2020):
        os.system("perl predict.pl .\\save_model\moel_1968_{} ..\\training_data\Du_Xiangjun\H3N2_seq_{} .\\save_ans\\ans_{} ".format(i,i+1,i+1))
    #perl modelGenerator.pl epitope_H3N2 aaIndex receptor_H3N2 ..\training_data\Du_Xiangjun\H3N2_seq_1968_2001 ..\training_data\Du_Xiangjun\seqCompare_H3N2_1968_2001 .\save_model\moel_1968_2001

if(choice =='5folds'):
    for i in range(5):
        os.system("perl predict.pl .\\save_model\moel_train_{} ..\\training_data_5folds\Du_Xiangjun\H3N2_seq_{} .\\save_ans\\ans_{} ".format(i,i,i))
