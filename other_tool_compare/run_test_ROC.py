from model_compare import *
import tensorflow as tf
import numpy as np
def main(feature_type,years1,years2):
    
    if(feature_type == 'Min-Shi Lee'):
        #(1) feature_type == 'Min-Shi Lee': don't need training
        for i in range(years1+1,years2+1):
            ans = Shi_Lee_test(i,i)
            a = np.array(ans)
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),a)
            for i in range(len(ans[0])):
                if(ans[0][i] != ans[1][i]):
                    print(i+2,end=",")
        print("\n")
    
    if(feature_type == 'Yu-Chieh Liao'):
        #(2) feature_type == 'Yu-Chieh Liao :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Liao_test_ROC(1968,i)
            a = np.array([ans_pre,ans_valid])
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),a)
        print("\n")

        
    if(feature_type == 'William Lees'):
        #(3) feature_type == 'William Lees :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = William_Lee_test_ROC(1968,i)
            a = np.array([ans_pre,ans_valid])
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),a)
        print("\n")
        
        
    if(feature_type == 'Peng Yousong'):
        #(4) feature_type == 'Peng Yousong :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Peng_test_ROC(1968,i)
            print(i)
            print(ans_pre)
            print(ans_valid)
            a = np.array([ans_pre,ans_valid])
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),a)
        print("\n")
        
    if(feature_type == 'Yuhua Yao'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Yao_test_ROC(1968,i)
            a = np.array([ans_pre,ans_valid])
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),a)
            print(ans_pre)
            print(ans_valid)
        print("\n")
        
        
    if(feature_type == 'IAV-CNN'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,loss,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,loss,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            ans_pre, ans_valid = IAV_CNN_test_ROC(1968,i)
            ans_valid = ans_valid.tolist()
            ans_valid = tf.keras.utils.to_categorical(ans_valid)

            ans_pre = ans_pre.tolist()
            print(type(ans_pre))
            print(type(ans_valid))
            np.save("./ROC/{}/{}.npy".format(feature_type,i+1),[ans_pre, ans_valid])
        file_T.close()
        file_V.close()
    
methods = ['Min-Shi Lee','Yu-Chieh Liao','William Lees','Peng Yousong','Yuhua Yao','IAV-CNN']
methods = ['IAV-CNN']
for method in methods:
    main(method,2001,2020)