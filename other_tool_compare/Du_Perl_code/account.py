import numpy as np
choice = '5folds' #'years'
if(choice == 'years'):
    for year in range(2001,2021):
        # origin H3N2 data
        csv = open("../training_data/Du_Xiangjun/seqCompare_H3N2_{}".format(year),"r")
        data_csv = []
        for line in csv.readlines():
            temp = line.split("\t")
            temp_data = [temp[0],temp[1],int(temp[2])]
            data_csv.append(temp_data)
        csv.close()
        #read the Du methods save_ans , count the data ['name1','name2',porbability]
        f = open("./save_ans/ans_{}".format(year),"r")
        data = []
        for line in f.readlines():
            temp = line.split("\t")
            temp_ans = 0
            if("M" in temp[2]):
                if(temp[2] == "Max"):
                    temp_ans = 1
                else:
                    temp_ans = 0
            else:
                temp_ans = float(temp[2])
            temp_data = [temp[0],temp[1],temp_ans]
            data.append(temp_data)
        f.close()
        
        #count the ans
        ans_pred = []
        ans_valid = []
        for i in data_csv:
            for k in data:
                if(i[0]==k[0] and i[1]==k[1] or i[0]==k[1] and i[1]==k[0]):
                    ans_pred.append(k[2])
                    ans_valid.append(i[2])
                    break
                else:
                    continue
        years_ans = np.array([ans_pred,ans_valid])
        np.save("../ROC/Du_Xiangjun/{}.npy".format(year),years_ans)
    
if(choice == '5folds'):
    for year in range(5):
        # origin H3N2 data
        csv = open("../training_data_5folds/Du_Xiangjun/seqCompare_H3N2_{}".format(year),"r")
        data_csv = []
        for line in csv.readlines():
            temp = line.split("\t")
            temp_data = [temp[0],temp[1],int(temp[2])]
            data_csv.append(temp_data)
        csv.close()
        #read the Du methods save_ans , count the data ['name1','name2',porbability]
        f = open("./save_ans/ans_{}".format(year),"r")
        data = []
        for line in f.readlines():
            temp = line.split("\t")
            temp_ans = 0
            if("M" in temp[2]):
                if(temp[2] == "Max"):
                    temp_ans = 1
                else:
                    temp_ans = 0
            else:
                temp_ans = float(temp[2])
            temp_data = [temp[0],temp[1],temp_ans]
            data.append(temp_data)
        f.close()
        
        #count the ans
        ans_pred = []
        ans_valid = []
        for i in data_csv:
            for k in data:
                if(i[0]==k[0] and i[1]==k[1] or i[0]==k[1] and i[1]==k[0]):
                    ans_pred.append(k[2])
                    ans_valid.append(i[2])
                    break
                else:
                    continue
        years_ans = np.array([ans_pred,ans_valid])
        np.save("../ROC/Du_Xiangjun/{}.npy".format(year),years_ans)