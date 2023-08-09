import numpy as np
import pandas as pd
import argparse
from model import *
import os





parser = argparse.ArgumentParser()
parser.add_argument('--number_filters', default=6, type=int)
parser.add_argument('--number_columns', default=14, type=int)
parser.add_argument('--filter_size', default=3, type=int, help='total number of training examples in a single batch')
parser.add_argument('--predict_data', required=True)
parser.add_argument('--model_path', required=True)
parser.add_argument('--outdir',required=True)
parser.add_argument('--type',required=True)
parser.add_argument('--seq_file', required=True, metavar='file')
args = parser.parse_args()



predict_data_0=args.predict_data
number_filters=args.number_filters
number_columns=args.number_columns
filter_size=args.filter_size
model_path=args.model_path
outdir_0=args.outdir+'/'
data_type=args.type

if os.path.exists(outdir_0):
    pass
else:
    os.makedirs(outdir_0)

if data_type=='H1N1':
    length=327
elif data_type=='H3N2':
    length=329

test_data=np.load(predict_data_0)


test_x=test_data[:,:,:14]



cnn = CNN(number_filters,filter_size,number_columns,length)



cnn.load_weights(model_path)
pred_y=cnn.predict(test_x)

pred_y_0=pred_y[:,1]
pred_y_0_np=np.array(pred_y_0)


np.save(outdir_0+'pred.npy',pred_y_0_np)
pred_y_0_np_01=np.where(pred_y_0_np >=0.5,1,0)

seq_all=pd.read_csv(args.seq_file)
seq_add=pd.read_csv(outdir_0+'pred_unselectedPairs_15var_5year.csv',sep='\t',index_col=0)
seq_add=seq_add.drop(['cha'], axis=1)
seq_all=seq_all[~seq_all.index.isin(seq_add.index)]
seq_all['predict']=pred_y_0_np_01
seq_all=pd.concat([seq_all,seq_add])

seq_all.to_csv(outdir_0+'pred.csv')
