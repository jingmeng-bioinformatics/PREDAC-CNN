import pandas as pd
filedata = pd.read_csv("./Antigen_H3N2_CNN_data.csv")
f = open("./junk.txt","w")
count = 0
len_data = filedata.shape[0]
#create indexlate data to save the csv , and count postive or negative sample
years_data = {}
years_data_count = {}
#range csv form begin to the end
for i in range(len_data):
    if(i%100==0):
        print("------{:.2f}%------".format(100*(1.0*i)/len_data),end="\r")
    pair1 = filedata['Strain1'][i]
    pair2 = filedata['Strain2'][i]
    if(pair1 == 'A/REASSORTANT/NYMC-X-161B'):
        pair1 = '2005'
    if(pair2 == 'A/REASSORTANT/NYMC-X-161B'):
        pair2 = '2005'
    pair1 = int(pair1[-4:])
    pair2 = int(pair2[-4:])
    index = i%5
    if(not index in years_data.keys()):
        years_data[index] = pd.DataFrame(columns=['Strain1','Strain2','Distance'])
        years_data[index] = years_data[index].append(filedata.iloc[i])
        years_data_count[index] = [0,0]
    else:
        years_data[index] = years_data[index].append(filedata.iloc[i])
1
        
#save all indexlate data
list_years = sorted(years_data.keys(),reverse=True)
ll = len(list_years)
for i in range(ll):  #   if i == ll-1 , data = all , we dont need it  ； data is i ,but name is i+1
    try:
        years_data[list_years[i]].to_csv("./data/{}.csv".format(list_years[i]),index=False)
    except:
        print("error")
            
f.close()
print("Congratulation!We finish it!")
    
