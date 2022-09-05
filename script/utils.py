import math
import numpy as np
import tensorflow as tf

H3N2_path = ".\\npy_del_single\\{}.npy"

H3N2_5folds_path_single = ".\\npy\\5folds\\{}.npy"
H3N2_5folds_path_double = ".\\npy\\5folds_double\\{}.npy"

def data_load_5folds_double(index):
    years_data = []  #index = [0,1,2,3,4]
    for i in index:
        try:
            years_data.append(np.load(H3N2_5folds_path_double.format(i)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data), axis=0)
    np.random.shuffle(years_data)
    x = years_data[:, :, :-1]
    x = x.reshape(-1, 329, 14, 1)
    y = years_data[:, 2, -1]
    y = tf.keras.utils.to_categorical(y)
    return x, y

def data_load_5folds_single(index):
    years_data = []  #index = [0,1,2,3,4]
    for i in index:
        try:
            years_data.append(np.load(H3N2_5folds_path_single.format(i)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data), axis=0)
    np.random.shuffle(years_data)
    x = years_data[:, :, :-1]
    x = x.reshape(-1, 329, 14, 1)
    y = years_data[:, 2, -1]
    y = tf.keras.utils.to_categorical(y)
    return x, y

def data_load_double_years(years_1, years_2):
    years_data = []
    for i in range(years_1, years_2 + 1):
        try:
            years_data.append(np.load(H3N2_path.format(i)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data), axis=0)
    np.random.shuffle(years_data)
    x = years_data[:, :, :-1]
    x = x.reshape(-1, 329, 14, 1)
    y = years_data[:, 2, -1]
    y = tf.keras.utils.to_categorical(y)
    return x, y

def data_load_test(years_1, years_2):
    years_data = []
    for i in range(years_1, years_2 + 1):
        try:
            years_data.append(np.load(H3N2_path.format(i)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data), axis=0)
    # np.random.shuffle(years_data)
    x = years_data[:, :, :-1]
    x = x.reshape(-1, 329, 14, 1)
    y = years_data[:, 2, -1]
    y = tf.keras.utils.to_categorical(y)
    return x, y

def data_load_single_years(years_1, years_2):
    #H3N2_path2 = "E:\\联想软件\\SCU\\0-CNN-Finall\\1-Data\\npy\\CNN_nodup_single\\{}.npy"
    H3N2_path2 = "E:\\联想软件\\SCU\\0-CNN-Finall\\1-Data\\del_test\\npy_del_single\\{}.npy"
    years_data = []
    for i in range(years_1, years_2 + 1):
        try:
            years_data.append(np.load(H3N2_path2.format(i)))
        except:
            continue
    years_data = np.concatenate(tuple(years_data), axis=0)
    # np.random.shuffle(years_data)
    x = years_data[:, :, :-1]
    x = x.reshape(-1, 329, 14, 1)
    y = years_data[:, 2, -1]
    y = tf.keras.utils.to_categorical(y)
    return x, y


def get_confusion_matrix(y_true, y_pred):
    """
    Calculates the confusion matrix from given labels and predictions.
    Expects tensors or numpy arrays of same shape.
    """
    TP, FP, TN, FN = 0, 0, 0, 0

    for i in range(y_true.shape[0]):
        if y_true[i] == 0 and y_pred[i] == 0:
            TN += 1
        elif y_true[i] == 0 and y_pred[i] == 1:
            FP += 1
        elif y_true[i] == 1 and y_pred[i] == 0:
            FN += 1
        elif y_true[i] == 1 and y_pred[i] == 1:
            TP += 1

    conf_matrix = [
        [TP, FP],
        [FN, TN]
    ]

    return conf_matrix


def get_accuracy(conf_matrix):
    """
    Calculates accuracy metric from the given confusion matrix.
    """
    TP, FP, FN, TN = conf_matrix[0][0], conf_matrix[0][1], conf_matrix[1][0], conf_matrix[1][1]
    return (TP + TN) / (TP + FP + FN + TN)


def get_precision(conf_matrix):
    """
    Calculates precision metric from the given confusion matrix.
    """
    TP, FP = conf_matrix[0][0], conf_matrix[0][1]

    if TP + FP > 0:
        return TP / (TP + FP)
    else:
        return 0


def get_recall(conf_matrix):
    """
    Calculates recall metric from the given confusion matrix.
    """
    TP, FN = conf_matrix[0][0], conf_matrix[1][0]

    if TP + FN > 0:
        return TP / (TP + FN)
    else:
        return 0


def get_f1score(conf_matrix):
    """
    Calculates f1-score metric from the given confusion matrix.
    """
    p = get_precision(conf_matrix)
    r = get_recall(conf_matrix)

    if p + r > 0:
        return 2 * p * r / (p + r)
    else:
        return 0


def get_mcc(conf_matrix):
    """
    Calculates Matthew's Correlation Coefficient metric from the given confusion matrix.
    """
    TP, FP, FN, TN = conf_matrix[0][0], conf_matrix[0][1], conf_matrix[1][0], conf_matrix[1][1]
    if TP + FP > 0 and TP + FN > 0 and TN + FP > 0 and TN + FN > 0:
        return (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    else:
        return 0


def evaluate(Y_pred, Y_real):
    conf_matrix = get_confusion_matrix(Y_real, Y_pred)
    precision = get_precision(conf_matrix)
    recall = get_recall(conf_matrix)
    fscore = get_f1score(conf_matrix)
    mcc = get_mcc(conf_matrix)
    val_acc = get_accuracy(conf_matrix)

    return precision, recall, fscore, mcc, val_acc


def Format_Y(Pre_y, Test_y):
    y1 = []
    y2 = []
    for i in range(Pre_y.shape[0]):
        if (Pre_y[i][0] >= 0.5):
            y1.append(1)
        else:
            y1.append(0)
        if (Test_y[i][0] == 1):
            y2.append(1)
        else:
            y2.append(0)
    return np.array(y1), np.array(y2)


def Format_Y2(Pre_y, Test_y):
    y1 = []
    y2 = []
    for i in range(Pre_y.shape[0]):
        if (Pre_y[i][0] >= 0.5):
            y1.append(0)
        else:
            y1.append(1)
        if (Test_y[i][0] == 1):
            y2.append(0)
        else:
            y2.append(1)
    return np.array(y1), np.array(y2)
