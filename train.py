import numpy as np
import pandas as pd
import argparse
from model import *
import os



parser = argparse.ArgumentParser()
parser.add_argument('--number_filters', default=6, type=int)
parser.add_argument('--number_columns', default=14, type=int)
parser.add_argument('--filter_size', default=3, type=int, help='total number of training examples in a single batch')
parser.add_argument('--epoch',default=100,type=int)
parser.add_argument('--train_data', required=True)
parser.add_argument('--test_data', required=True)
parser.add_argument('--outdir',required=True)
parser.add_argument('--type',required=True)
args = parser.parse_args()


train_data_0=args.train_data
test_data_0=args.test_data
number_filters=args.number_filters
number_columns=args.number_columns
filter_size=args.filter_size
epoch=args.epoch
data_type=args.type
outdir_0=args.outdir+'/'


if os.path.exists(outdir_0):
    pass
else:
    os.makedirs(outdir_0)

if data_type=='H1N1':
    length=327
elif data_type=='H3N2':
    length=329

train_data=np.load(train_data_0)
test_data=np.load(test_data_0)

train_x_0=train_data[:,:,:14]
train_x_copy=np.roll(train_x_0,-7,axis=2)
train_x=np.concatenate((train_x_0,train_x_copy),axis=0)

train_y_0_0=train_data[:,0,14]

train_y_0=np.concatenate((train_y_0_0,train_y_0_0))
test_x=test_data[:,:,:14]
test_y_0=test_data[:,0,14]
train_y=np.zeros((train_x.shape[0],2))
test_y=np.zeros((test_x.shape[0],2))
for i in range(train_x.shape[0]):
    train_y[i,0]=1-train_y_0[i]
    train_y[i,1]=train_y_0[i]

for i in range(test_x.shape[0]):
    test_y[i,0]=1-test_y_0[i]
    test_y[i,1]=test_y_0[i]


checkpoiner=tf.keras.callbacks.ModelCheckpoint(filepath=outdir_0+'/model_{epoch:02d}',monitor='val_loss',save_weights_only=True,verbose=1)
cnn = CNN(number_filters,filter_size,number_columns,length)
history_callback=cnn.fit(
        x=train_x, y=train_y,
        batch_size=128, 
        epochs=epoch, 
        verbose=1, 
        callbacks=[checkpoiner],
        shuffle=True,
        validation_data=(test_x,test_y))

TP_all=[]
FP_all=[]
FN_all=[]
TN_all=[]

for i in range(epoch):
    i+=1
    i='%02d' %i
    cnn.load_weights(outdir_0+'/model_'+str(i))
    pred_y=cnn.predict(test_x)
    TP, FP, TN, FN=get_confusion_matrix(test_y, pred_y)
    TP_all.append(TP)
    FP_all.append(FP)
    TN_all.append(TN)
    FN_all.append(FN)


new_dict=history_callback.history.copy()
new_dict['TN']=TN_all
new_dict['FP']=FP_all
new_dict['FN']=FN_all
new_dict['TP']=TP_all



pandas.DataFrame(new_dict).to_csv(outdir_0 + '/log.csv')


