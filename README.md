# PREDAC-CNN

PREDAC-CNN is based on the CNN to predict the antigenic variants of seasonal influenza A viruses. PREDAC-CNN constructs a spatially oriented representation of the HA1 sequence adapted for the convolutional architecture, which can explore the interactions of the amino acid sites in the context sequence. Moreover, rather than the redundant amino acid embeddings, PREDAV-CNN takes into account only physicochemical features determining antigenicity of influenza viruses. Together, PREDAC-CNN can effectively capture the dependencies in the context sequence, and investigate the combinatorial contributions of point mutations in the HA protein to the antigenic variants.

## Prerequisites

    tensorflow==2.4.0
    python==3.7
    pandas==1.3.5


## Getting started
1 Please prepare the input file. It is required that the input file contains 'name1','name2','seq1','seq2','year1','year2' and 'label' for training data and testing data or 'name1','name2','seq1','seq2','year1','year2' for prediction.<br>

Example of the training data and testing data:



![image](https://github.com/jingmeng-bioinformatics/PREDAC-CNN/assets/35085665/ce68d776-0876-4ad0-8866-2d50954b941c)



Example of the input data for prediction:



![image](https://github.com/jingmeng-bioinformatics/PREDAC-CNN/assets/35085665/42abde09-d86e-4cb9-8840-592327329e33)




2 Run matrix_generate.py to generate the input matrix from the input file (train_data,test_data or predict_data):

    python matrix_generate.py
    --aaindex_file aaindex_feature_H1N1.txt (or aaindex_feature_H3N2.txt)
    --seq_file /path/input_file
    --type training (or testing or predict)
    --dir /path/save_dir/


3 Run train.py to train a CNN model:

    python train.py
    --train_data /path/save_dir/train_data
    --test_data /path/save_dir/test_data
    --outdir /path/model_dir/
    --type H1N1 (or H3N2)

4 Run predict.py to predict the antigenic variants of the two sequences:

    python predict.py
    --predict_data /path/save_dir/predict_data
    --seq_file /path/input_file
    --model_path /path/model_dir/model
    --outdir /path/predict_dir/
    --type H1N1 (or H3N2)

Example of the output file:

![image](https://github.com/jingmeng-bioinformatics/PREDAC-CNN/assets/35085665/29ecd1b6-3c93-492e-8503-cd357c8706b8)

5 Run mcl_cluster.py to predict the antigenic variants of the two sequences:

    python3 mcl_cluster.py
    --file /path/predict_dir/predict_data
    --dir /path/save_dir/
    --type H1N1 (or H3N2)




Jing Meng<br>

jing.mengrabbit@outlook.com<br>
