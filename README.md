# PREDAV-CNN

PREDAV-CNN is based on the CNN to predict the antigenic variants of seasonal influenza A viruses. PREDAV-CNN constructs a spatially oriented representation of the HA1 sequence adapted for the convolutional architecture, which can explore the interactions of the amino acid sites in the context sequence. Moreover, rather than the redundant amino acid embeddings, PREDAV-CNN takes into account only physicochemical features determining antigenicity of influenza viruses. Together, PREDAV-CNN can effectively capture the dependencies in the context sequence, and investigate the combinatorial contributions of point mutations in the HA protein to the antigenic variants.

## Prerequisites

    tensorflow==2.4.0
    python==3.7
    pandas==1.3.5


## Getting started
1 Please prepare the input file. It is required that the input file contains 'seq1','seq2' and 'label' for training data and testing data or 'seq1'and 'seq2' for prediction.<br>

Example of the training data and testing data:

    label    seq_1            seq_2
    1        QKLPGNDNST...    QKLPGNDNSS...
    0        -KLPGNDNSS...    -KLPGNDNT...
    1        QKLPGIDNSN...    QKLPGIDNSS...
    0        QKLPGNDNTS...    QKLPGNDNSS...

Example of the input data for prediction:

    seq_1	        seq_2
    QKLPGNDNST...	QKLPGNDNSS...
    -KLPGNDNS...	-KLPGNDNT...
    QKLPGIDNSN...	QKLPGIDNSS...
    QKLPGNDNTS...	QKLPGNDNSS...


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

4 Run predict.py to the antigenic variants of the two sequences:

    python predict.py
    --predict_data /path/save_dir/predict_data
    --seq_file /path/input_file
    --model_path /path/model_dir/model
    --outdir /path/predict_dir/
    --type H1N1 (or H3N2)

Example of the output file:

    seq_1	        seq_2	        predict
    QKLPGNDNST...	QKLPGNDNSS...	0
    -KLPGNDNS...	-KLPGNDNT...	1
    QKLPGIDNSN...	QKLPGIDNSS...	0
    QKLPGNDNTS...	QKLPGNDNSS...	1



Jing Meng<br>

jing.mengrabbit@outlook.com<br>
