from model_compare import *

def main(feature_type,years1,years2):
    if(feature_type == 'Min-Shi Lee'):
        #(1) feature_type == 'Min-Shi Lee': don't need training
        file = open("./log/{}.csv".format(feature_type),'w')
        file.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            ans = Shi_Lee_method(i,i)
            print(type(ans))
            file.write(ans)
        file.close()
    
    if(feature_type == 'Yu-Chieh Liao'):
        #(2) feature_type == 'Yu-Chieh Liao :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            T_ans, V_ans = Liao_method(1968,i)
            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'William Lees'):
        #(3) feature_type == 'William Lees :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            T_ans, V_ans = William_Lee_method(1968,i)
            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'Peng Yousong'):
        #(4) feature_type == 'Peng Yousong :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            T_ans, V_ans = Peng_method(1968,i)
            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'Yuhua Yao'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            T_ans, V_ans = Yao_method(1968,i)
            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
        
        
    if(feature_type == 'IAV_CNN_method'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,loss,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,loss,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            T_ans, V_ans = IAV_CNN_method(1968,i)
            #      print('V_loss %.3f\tV_acc %.3f\tV_pre %.3f\tV_rec %.3f\tV_fscore %.3f\tV_mcc %.3f'%(val_loss, val_acc, precision, recall, fscore, mcc))

            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
    
methods = ['Min-Shi Lee','Yu-Chieh Liao','William Lees','Peng Yousong','Yuhua Yao']#'IAV_CNN_method'
methods = ['IAV_CNN_method']
for method in methods:
    try:
        main(method,2000,2021)
    except:
        continue