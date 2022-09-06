from keras.layers import Conv2D, MaxPool2D, Dense, Flatten, Dropout, LayerNormalization
from keras.models import Sequential
import numpy as np
from keras.layers.core import Activation
from tensorflow.keras.optimizers import SGD, Adam, RMSprop
from keras.callbacks import ReduceLROnPlateau
from utils import *

def CNN(lr=1e-3):
    lenet = Sequential()
    # 表示我们的网络将学习6个滤波器 每个滤波器的大小都是3×3，步长为1 padding='',
    lenet.add(Conv2D(6, kernel_size=(3, 14), strides=1, input_shape=(329, 14, 1)))  # 329 --> 166
    lenet.add(Activation("relu"))
    lenet.add(MaxPool2D(pool_size=(4, 1), strides=2))

    lenet.add(Conv2D(12, kernel_size=(7, 1), strides=3, input_shape=(162, 1, 6)))
    lenet.add(Activation("relu"))
    lenet.add(MaxPool2D(pool_size=(2, 1), strides=2))
    lenet.add(Flatten())
    lenet.add(Dense(120))
    lenet.add(Activation("relu"))
    lenet.add(Dropout(0.25))
    # 84个神经元 全连接网络
    lenet.add(Dense(84))
    lenet.add(Activation("relu"))
    lenet.add(Dense(2, activation='softmax'))  # 10个类别的softmax分类器
    lenet.compile(optimizer=Adam(lr=lr, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0, clipvalue=0.5),
                  loss='categorical_crossentropy', metrics=['accuracy'])
    return lenet


def run_years(years1, years2):
    file_orig = open("CNN_years.csv", "w")
    file_orig.write("years,epoch,precision,recall,fscore,mcc,val_acc\n")
    file_orig.close()
    for i in range(years1, years2):
        if (True):
            best_acc = 0
            file_orig = open("CNN_years.csv", "a")
            Train_x, Train_y = data_load_double_years(1968, i)
            Test_x, Test_y = data_load_single_years(i + 1, i + 1)
            cnn = CNN(lr=0.001)
            print(cnn.summary())
            reduce_Ir = ReduceLROnPlateau(monitor='loss', patience=3, mode='auto', min_lr=0.0001, factor=0.5)
            for k in range(20):
                cnn.fit(Train_x, Train_y, callbacks=[reduce_Ir], batch_size=32, epochs=4)
                Pre_y = cnn.predict(Test_x)
                temp = Format_Y2(Pre_y, Test_y)
                precision, recall, fscore, mcc, val_acc = evaluate(temp[0], temp[1])
                score1 = cnn.evaluate(Test_x, Test_y)
                if (val_acc > best_acc):
                    best_acc = val_acc
                    cnn.save_weights("./model-onlysingle/{}.h5".format(i + 1))
                    best = "{},".format(i + 1) + "{},".format(k) + "{},{},{},{},{},\n".format(precision, recall, fscore,
                                                                                              mcc, val_acc)
            file_orig.write(best)
            file_orig.close()

        else:
            continue

def run_5folds():
    file_orig = open("CNN_5folds.csv", "w")
    file_orig.write("years,epoch,precision,recall,fscore,mcc,val_acc\n")
    file_orig.close()
    for i in range(1):
        index_train = []
        index_test = []
        for j in range(5): #choose i as the test data
            if(i==j):
                index_test.append(j)
            else:
                index_train.append(j) #others as the train data

        if (True):
            best_acc = 0
            file_orig = open("CNN_5folds.csv", "a")
            Train_x, Train_y = data_load_5folds_double(index_train)
            Test_x, Test_y = data_load_5folds_single(index_test)
            cnn = CNN()
            print(cnn.summary())
            reduce_Ir = ReduceLROnPlateau(monitor='loss', patience=3, mode='auto', min_lr=0.0001, factor=0.5)
            for k in range(40):
                cnn.fit(Train_x, Train_y, callbacks=[reduce_Ir], batch_size=32, epochs=4)
                Pre_y = cnn.predict(Test_x)
                temp = Format_Y2(Pre_y, Test_y)
                precision, recall, fscore, mcc, val_acc = evaluate(temp[0], temp[1])
                score1 = cnn.evaluate(Test_x, Test_y)
                if (val_acc > best_acc):
                    best_acc = val_acc
                    cnn.save_weights("./model_5folds/{}.h5".format(i + 1))
                    best = "{},".format(i + 1) + "{},".format(k) + "{},{},{},{},{},\n".format(precision, recall, fscore,
                                                                                              mcc, val_acc)
            file_orig.write(best)
            file_orig.close()

        else:
            continue

def run_test_years(years):
    f = open("./Finall_ans.csv", "a")
    f2 = open("./years_ROC.data","a")
    single_x, single_y = data_load_single_years(years, years)
    cnn = CNN()
    cnn.load_weights("./model-onlysingle/{}.h5".format(years))
    pre_sy = cnn.predict(single_x)
    temp2 = Format_Y2(pre_sy, single_y)
    precision, recall, fscore, mcc, val_acc = evaluate(temp2[0], temp2[1])
    y_pre, y_test = Format_Y2(pre_sy, single_y)

    a = np.array([pre_sy, single_y])
    np.save("./ROC/H3N2-CNN_only_single/{}.npy".format(years),a)
    f2.write("{}\n".format(str(y_pre)))
    f2.write("{}\n".format(str(y_test)))
    print(precision, recall, fscore, mcc, val_acc)
    print(y_pre)
    print(y_test)
    for i in range(len(y_pre)):
        if(y_pre[i] != y_test[i]):
            print(i+2,end=",")
    f.write("{},{},{},{},{},{}\n".format(years,precision,recall,fscore,mcc,val_acc))
    f2.close()
    f.close()

def main():
    #run_years(2000, 2020)
    for i in range(2001,2021):
        run_test_years(i)

def run_test():
    f = open("./Finall_ans.csv", "w")
    f.write("Years,precision,recall,fscore,mcc,val_acc\n")
    f.close()
    for i in range(2001, 2021):
        run_test_years(years=i)
        
#run_test()
cnn = CNN()
print(cnn.summary())
