from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Conv1D
import tensorflow as tf
from tensorflow.keras.callbacks import ReduceLROnPlateau
import math
import numpy as np
import os
from tensorflow.keras import backend as K

def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true[:,1] * y_pred[:,1], 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true[:,1], 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true[:,1] * y_pred[:,1], 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred[:,1], 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

def get_confusion_matrix(y_true, y_pred):
    """
    Calculates the confusion matrix from given labels and predictions.
    Expects tensors or numpy arrays of same shape.
    """
    TP, FP, TN, FN = 0, 0, 0, 0
    for i in range(y_pred.shape[0]):
        if y_pred[i,1] >= 0.5:
            y_pred[i,1]=1
            y_pred[i,0]=0
        else:
            y_pred[i,1]=0
            y_pred[i,0]=1

    TP=sum(y_true[:,1] * y_pred[:,1])
    TN=sum(y_true[:,0] * y_pred[:,0])
    FP=sum(y_true[:,0] * y_pred[:,1])
    FN=sum(y_true[:,1] * y_pred[:,0])


    return TP, FP, TN, FN


def CNN(number_filters,filter_size,number_columns,length,lr=1e-3):
    model = tf.keras.Sequential([
    tf.keras.layers.Conv1D(number_filters, kernel_size=filter_size, padding='same',strides=1, input_shape=(length, number_columns),activation='relu'),
    tf.keras.layers.Conv1D(number_filters, kernel_size=filter_size, padding='same',strides=1, input_shape=(length, number_filters),activation='relu'),
    tf.keras.layers.MaxPooling1D(2,2),
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(128,activation='relu'),
    tf.keras.layers.Dropout(0.25),
    tf.keras.layers.Dense(64,activation='relu'),
    tf.keras.layers.Dense(2, activation='softmax')])
    model.compile(optimizer=tf.keras.optimizers.Adam(lr=lr, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0, clipvalue=0.5),
                  loss='categorical_crossentropy', metrics=['accuracy',f1_m])
    return model