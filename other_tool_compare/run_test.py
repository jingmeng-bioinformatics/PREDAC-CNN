from model_compare import *

def main(feature_type,years1,years2):
    
    if(feature_type == 'Min-Shi Lee'):
        #(1) feature_type == 'Min-Shi Lee': don't need training
        for i in range(years1,years2):
            ans = Shi_Lee_test(i,i)
            print(ans[0])
            print(len(ans[0]))
            print(ans[1])
            print(len(ans[1]))
            for i in range(len(ans[0])):
                if(ans[0][i] != ans[1][i]):
                    print(i+2,end=",")
        print("\n")
    
    if(feature_type == 'Yu-Chieh Liao'):
        #(2) feature_type == 'Yu-Chieh Liao :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Liao_test(1968,i)
            print(ans_pre)
            print(ans_valid)
            
            for i in range(len(ans_pre)):
                if(ans_pre[i] != ans_valid[i]):
                    print(i+2,end=",")

        print("\n")

        
    if(feature_type == 'William Lees'):
        #(3) feature_type == 'William Lees :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = William_Lee_test(1968,i)
            print(ans_pre)
            print(ans_valid)
            for i in range(len(ans_pre)):
                if(ans_pre[i] != ans_valid[i]):
                    print(i+2,end=",")
        print("\n")
        
        
    if(feature_type == 'Peng Yousong'):
        #(4) feature_type == 'Peng Yousong :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Peng_test(1968,i)
            print(ans_pre)
            print(ans_valid)
            for i in range(len(ans_pre)):
                if(ans_pre[i] != ans_valid[i]):
                    print(i+2,end=",")
        print("\n")
        
    if(feature_type == 'Yuhua Yao'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        for i in range(years1,years2):
            ans_pre, ans_valid = Yao_test(1968,i)
            print(ans_pre)
            print(ans_valid)
            for i in range(len(ans_pre)):
                if(ans_pre[i] != ans_valid[i]):
                    print(i+2,end=",")
        print("\n")
        
        
    if(feature_type == 'IAV_CNN_method'):
        #(5) feature_type == 'Yuhua Yao :': training is years1~years2 ,test is years2+1
        file_T = open("./log/{}_train.csv".format(feature_type),'w')
        file_T.write("Years,loss,T_acc,T_pre,T_rec,T_fscore,T_mcc\n")
        
        file_V = open("./log/{}_valid.csv".format(feature_type),'w')
        file_V.write("Years,loss,V_acc,V_pre,V_rec,V_fscore,V_mcc\n")
        for i in range(years1,years2):
            ans_pre, ans_valid = IAV_CNN_method(1968,i)
            #      print('V_loss %.3f\tV_acc %.3f\tV_pre %.3f\tV_rec %.3f\tV_fscore %.3f\tV_mcc %.3f'%(val_loss, val_acc, precision, recall, fscore, mcc))

            file_T.write("{}-{},".format(1968,i)+T_ans)
            file_V.write("{},".format(i+1)+V_ans)
        file_T.close()
        file_V.close()
    
methods = ['Min-Shi Lee','Yu-Chieh Liao','William Lees','Peng Yousong','Yuhua Yao']
methods = ['Yu-Chieh Liao','William Lees','Peng Yousong']
for method in methods:
    try:
        main(method,2001,2020)
    except:
        continue