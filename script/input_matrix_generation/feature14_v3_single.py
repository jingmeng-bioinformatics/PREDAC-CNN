import pandas as pd
import numpy as np


aj_dic = {
    'A': {'access': 0.09, 'charge': 0.0, 'hydro': 0.23, 'hyindex': 0, 'polar': 0,'volume':0.47, 'num': 1}, 
    'R': {'access': 0.91, 'charge': 1, 'hydro': 0.23, 'hyindex': 1, 'polar': 1.0,'volume':0.60, 'num': 2}, 
    'N': {'access': 0.57, 'charge': 0.0, 'hydro': 0.02, 'hyindex': 0.5, 'polar': 0.07,'volume':0.13, 'num': 3}, 
    'D': {'access': 0.41, 'charge': -1, 'hydro': 0.17, 'hyindex': 0.25, 'polar': 0.96,'volume':0.00, 'num': 4}, 
    'C': {'access': 0.01, 'charge': 0.0, 'hydro': 0.40, 'hyindex': 0, 'polar': 0.03,'volume':0.23, 'num': 5},  
    'Q': {'access': 0.67, 'charge': 0.0, 'hydro': 0.00, 'hyindex': 0.5, 'polar': 0.07,'volume':0.32, 'num': 6}, 
    'E': {'access': 0.39, 'charge': -1, 'hydro': 0.18, 'hyindex': 0.25, 'polar': 0.96,'volume':0.26, 'num': 7},
    'G': {'access': 0.06, 'charge': 0.0, 'hydro': 0.03, 'hyindex': 0, 'polar': 0.0,'volume':0.12, 'num': 8}, 
    'H': {'access': 0.32, 'charge': 0.0, 'hydro': 0.23, 'hyindex': 0.25, 'polar': 0.99,'volume':0.31, 'num': 9}, 
    'I': {'access': 0, 'charge': 0.0, 'hydro': 0.84, 'hyindex': 0, 'polar': 0,'volume':0.88, 'num': 10}, 
    'L': {'access': 0.06, 'charge': 0.0, 'hydro': 0.58, 'hyindex': 0.0, 'polar': 0,'volume':1.00, 'num': 11}, 
    'K': {'access': 1, 'charge': 1, 'hydro': 0.43, 'hyindex': 0.5, 'polar': 0.95,'volume':0.74, 'num': 12}, 
    'M': {'access': 0.16, 'charge': 0.0, 'hydro': 0.45, 'hyindex': 0, 'polar': 0.03,'volume':0.53, 'num': 13}, 
    'F': {'access': 0.08, 'charge': 0.0, 'hydro': 0.76, 'hyindex': 0.0, 'polar': 0.01,'volume':0.70, 'num': 14}, 
    'P': {'access': 0.41, 'charge': 0.0, 'hydro': 0.74, 'hyindex': 0, 'polar': 0.03,'volume':0.61, 'num': 15}, 
    'S': {'access': 0.33, 'charge': 0.0, 'hydro': 0.02, 'hyindex': 0.25, 'polar': 0.03,'volume':0.13, 'num': 16}, 
    'T': {'access': 0.37, 'charge': 0.0, 'hydro': 0.02, 'hyindex': 0.25, 'polar': 0.03,'volume':0.34, 'num': 17}, 
    'W': {'access': 0.18, 'charge': 0.0, 'hydro': 1.00, 'hyindex': 0.25, 'polar': 0.04,'volume':0.65, 'num': 18}, 
    'Y': {'access': 0.53, 'charge': 0.0, 'hydro': 0.71, 'hyindex': 0.25, 'polar': 0.03,'volume':0.65, 'num': 19}, 
    'V': {'access': 0, 'charge': 0.0, 'hydro': 0.50, 'hyindex': 0.00, 'polar': 0,'volume':0.77, 'num': 20}, 
    }

mark_arr = [50,57,121,122,124,129,131,132,133,135,137,140,142,143,144,145,146,152,155,156,157,157,158,159,160,164,172,173,188,189,190,193,196,197,207,208,216,217,219,225,226,240,244,260,275,276,278,279]
class str_to_num():
    def __init__(self):
        self.seq1 = []
        self.seq2 = []
        self.seq_class = []
        self.error_code = -1
        #version2
        self.all_data = []


    
    def read_from_csv(self,filename):
        print("starting read data....")
        csv_data = pd.read_csv(filename+".csv", encoding= 'utf-8')
        self.seq1 = csv_data['Seq1']
        self.seq2 = csv_data['Seq2']
        self.seq_class = csv_data['Class']

        assert(len(self.seq1) == len(self.seq2))
        assert(len(self.seq1) == len(self.seq_class))

        print("read data done")
    

    #version 2
    def generate_seq(self, seq1, seq2 ,clas):
        #this is a two-dimensional array
        one_seq = []
        if seq1 and seq2:
            for i in range(len(seq1)):
                #this is one-dimensional array, to store a character info
                one_char = [-1]*15
                if seq1[i] in aj_dic:
                    one_char[0] = aj_dic[seq1[i]]['num'] / 20
                    one_char[1] = aj_dic[seq1[i]]['access']
                    one_char[2] = aj_dic[seq1[i]]['charge']
                    one_char[3] = aj_dic[seq1[i]]['hydro']
                    one_char[4] = aj_dic[seq1[i]]['hyindex']
                    one_char[5] = aj_dic[seq1[i]]['polar']
                    one_char[6] = aj_dic[seq1[i]]['volume']
                    
                if seq2[i] in aj_dic:
                    one_char[7] = aj_dic[seq2[i]]['num'] / 20
                    one_char[8] = aj_dic[seq2[i]]['access']
                    one_char[9] = aj_dic[seq2[i]]['charge']
                    one_char[10] = aj_dic[seq2[i]]['hydro']
                    one_char[11] = aj_dic[seq2[i]]['hyindex']
                    one_char[12] = aj_dic[seq2[i]]['polar']
                    one_char[13] = aj_dic[seq2[i]]['volume']
                one_char[14] = clas
                one_seq.append(one_char)    
        else:
            print("sequence is null...")
        return one_seq
        # assert(len(one_seq) == 325)
    #数组存顺序：seq1, seq2, hd,ha,pc,nc,h,mark
    def do_generate(self,filename):
        #self.read_from_csv("./CNN_nodup_single/"+filename)
        self.read_from_csv("./CNN_nodup_double/"+filename)
        print('start generating seq....')
        for i in range(len(self.seq1)):
            if i % 1000 == 0:
                print("i = ", i)
                #print(self.all_data)
            self.all_data.append(self.generate_seq(self.seq1[i], self.seq2[i],self.seq_class[i]))
        print('generate seq done')

    def save_to_npy(self,filename):
        arr = np.array(self.all_data)
        #np.save('./npy_del_single/'+filename+'.npy', arr)
        np.save('./npy_del_double/'+filename+'.npy', arr)


if __name__ == "__main__":
    # print("hello world")
    for i in range(1968,2022):
        try:
            s = str_to_num()
            s.do_generate(str(int(i)))
            # s.save_file()
            s.save_to_npy(str(int(i)))
        except:
            continue
        
