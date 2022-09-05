from model_compare import *

def main(feature_type):
    
    if(feature_type == 'Min-Shi Lee'):
        #(1) feature_type == 'Min-Shi Lee': don't need training
        file = open("./log_5folds/{}.csv".format(feature_type),'w')
        file.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(5):
            valid = []
            valid.append(i)
            ans = Shi_Lee_5folds_method(valid)
            file.write(ans)
        file.close()
    
    if(feature_type == 'Yu-Chieh Liao'):
        #(2) feature_type == 'Yu-Chieh Liao :': training is [0,1,2,3] ,test is [4]
        file_T = open("./log_5folds/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log_5folds/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        
        for i in range(5):
            train = []
            valid = []
            for k in range(5):
                if(k!=i):
                    train.append(k)
                else:
                    valid.append(i)
            T_ans, V_ans = Liao_5folds_method(train,valid)
            file_T.write("{},".format(str(train))+T_ans)
            file_V.write("{},".format(valid)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'William Lees'):
        #(3) feature_type == 'William Lees :': training is [0,1,2,3] ,test is [4
        file_T = open("./log_5folds/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log_5folds/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")

        for i in range(5):
            train = []
            valid = []
            for k in range(5):
                if(k!=i):
                    train.append(k)
                else:
                    valid.append(i)
            T_ans, V_ans = William_Lee_5folds_method(train,valid)
            file_T.write("{},".format(str(train))+T_ans)
            file_V.write("{},".format(valid)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'Peng Yousong'):
        #(4) feature_type == 'Peng Yousong :': training is years1~years2 ,test is years2+1
        file_T = open("./log_5folds/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        file_V = open("./log_5folds/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")

        for i in range(5):
            train = []
            valid = []
            for k in range(5):
                if(k!=i):
                    train.append(k)
                else:
                    valid.append(i)
            T_ans, V_ans = Peng_5folds_method(train,valid)
            file_T.write("{},".format(str(train))+T_ans)
            file_V.write("{},".format(valid)+V_ans)
        file_T.close()
        file_V.close()
        
    if(feature_type == 'Yuhua Yao'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log_5folds/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        file_V = open("./log_5folds/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")

        for i in range(5):
            train = []
            valid = []
            for k in range(5):
                if(k!=i):
                    train.append(k)
                else:
                    valid.append(i)
            T_ans, V_ans = Yao_5folds_method(train,valid)
            file_T.write("{},".format(str(train))+T_ans)
            file_V.write("{},".format(valid)+V_ans)
        file_T.close()
        file_V.close()
        
        
    if(feature_type == 'IAV_CNN_5folds_method'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log_5folds/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,loss,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log_5folds/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,loss,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(5):
            train = []
            valid = []
            for k in range(5):
                if(k!=i):
                    train.append(k)
                else:
                    valid.append(i)
            T_ans, V_ans = IAV_CNN_5folds_method(train,valid)
            file_T.write("{},".format(str(train))+T_ans)
            file_V.write("{},".format(valid)+V_ans)
        file_T.close()
        file_V.close()
#'Min-Shi Lee',
methods = ['Yu-Chieh Liao','William Lees','Peng Yousong','Yuhua Yao','IAV_CNN_5folds_method']
methods = ['Yu-Chieh Liao']
for method in methods:
    main(method)