import pandas as pd
import numpy as np
import argparse





parser = argparse.ArgumentParser()

parser.add_argument('--aaindex_file', required=True, metavar='file')
parser.add_argument('--seq_file', required=True, metavar='file')
parser.add_argument('--type', required=True,metavar='predict or training or testing')
parser.add_argument('--dir', required=True, metavar='file')

args = parser.parse_args()



class str_to_num():
    def __init__(self):
        self.seq1 = []
        self.seq2 = []
        #version2
        self.all_data = []

    def read_from_csv(self,csv_data,type='predict'):
        print("starting read data....")
        self.seq1 = csv_data['seq_1']
        self.seq2 = csv_data['seq_2']
        if type=='predict':
            pass
        else:
            self.label = csv_data['label']
            assert(len(self.seq1) == len(self.label))
        assert(len(self.seq1) == len(self.seq2))


        print("read data done")

    def generate_seq_test(self,seq1, seq2 ):
            #this is a two-dimensional array
            one_char = np.zeros((len(seq1),len(seq1[0]),14))
            for i in range(len(seq1)):
                for j in range(len(seq1[i])):
                    one_char[i][j][0] = aj_dic[seq1[i][j]]['number']
                    one_char[i][j][1] = aj_dic[seq1[i][j]]['access']
                    one_char[i][j][2] = aj_dic[seq1[i][j]]['charge']
                    one_char[i][j][3] = aj_dic[seq1[i][j]]['hydro']
                    one_char[i][j][4] = aj_dic[seq1[i][j]]['hyindex']
                    one_char[i][j][5] = aj_dic[seq1[i][j]]['polar']
                    one_char[i][j][6] = aj_dic[seq1[i][j]]['volume']
                    one_char[i][j][7] = aj_dic[seq2[i][j]]['number']
                    one_char[i][j][8] = aj_dic[seq2[i][j]]['access']
                    one_char[i][j][9] = aj_dic[seq2[i][j]]['charge']
                    one_char[i][j][10] = aj_dic[seq2[i][j]]['hydro']
                    one_char[i][j][11] = aj_dic[seq2[i][j]]['hyindex']
                    one_char[i][j][12] = aj_dic[seq2[i][j]]['polar']
                    one_char[i][j][13] = aj_dic[seq2[i][j]]['volume']
            return one_char

    def generate_seq(self,seq1, seq2 ,label):
            #this is a two-dimensional array
            one_char = np.zeros((len(seq1),len(seq1[0]),15))
            for i in range(len(seq1)):
                for j in range(len(seq1[i])):
                    one_char[i][j][0] = aj_dic[seq1[i][j]]['number']
                    one_char[i][j][1] = aj_dic[seq1[i][j]]['access']
                    one_char[i][j][2] = aj_dic[seq1[i][j]]['charge']
                    one_char[i][j][3] = aj_dic[seq1[i][j]]['hydro']
                    one_char[i][j][4] = aj_dic[seq1[i][j]]['hyindex']
                    one_char[i][j][5] = aj_dic[seq1[i][j]]['polar']
                    one_char[i][j][6] = aj_dic[seq1[i][j]]['volume']
                    one_char[i][j][7] = aj_dic[seq2[i][j]]['number']
                    one_char[i][j][8] = aj_dic[seq2[i][j]]['access']
                    one_char[i][j][9] = aj_dic[seq2[i][j]]['charge']
                    one_char[i][j][10] = aj_dic[seq2[i][j]]['hydro']
                    one_char[i][j][11] = aj_dic[seq2[i][j]]['hyindex']
                    one_char[i][j][12] = aj_dic[seq2[i][j]]['polar']
                    one_char[i][j][13] = aj_dic[seq2[i][j]]['volume']
                    one_char[i][j][14] = label[i]
            return one_char

    def do_generate(self,type='predict'):
        print('start generating seq....')
        if type=='predict':
            self.all_data=self.generate_seq_test(self.seq1, self.seq2)
        else:
            self.all_data=self.generate_seq(self.seq1, self.seq2,self.label)
        print('generate seq done')


    def save_to_npy(self,dir,filename):
        arr = np.array(self.all_data)
        np.save(dir+'/'+filename+'.npy', arr)


##################################################


seq_all=pd.read_csv(args.seq_file,sep='\t',index_col=0)
aaindex=pd.read_csv(args.aaindex_file,sep='\t',index_col=0)
dir=args.dir
all_dict=aaindex.T.to_dict()
aj_dic=all_dict
type=args.type
##################################################

s = str_to_num()
s.read_from_csv(seq_all,type)
s.do_generate(type)
s.save_to_npy(dir,type)