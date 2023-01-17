# PREDAV-CNN

CNNantigentic is based on the CNN to predict the antigenic relationship of influenza A/H3N2 virus. CNNantigenic constructs a spatially oriented representation of the HA1 sequence adapted for the convolutional architecture, which can explore the interactions of the amino acids in the context sequence. Moreover, rather than the redundant amino acid embeddings, CNNantigenic takes into account only physicochemical features determining antigenicity of influenza A/H3N2 virus.Together, CNNantigenic can effectively extract the features in the context sequence from local to global views, and investigate the combinatorial contributions of point mutations in the HA protein to the antigenicity.

## Prerequisites

tensorflow==2.4.0
python==3.7
pandas==1.3.5


## Getting started
1 Please prepare the input file. It is required that input file has to contain 'seq1','seq2' and 'label' for training data and testing data or 'seq1'and 'seq2' for predict.
Example of the training data and testing data

label	seq_1	seq_2
1	QKLPGNDNST...	QKLPGNDNSS...
0	-KLPGNDNS...	-KLPGNDNT...
1	QKLPGIDNSN...	QKLPGIDNSS...
0	QKLPGNDNTS...	QKLPGNDNSS...

Example of the predict data

seq_1	seq_2
QKLPGNDNST...	QKLPGNDNSS...
-KLPGNDNS...	-KLPGNDNT...
QKLPGIDNSN...	QKLPGIDNSS...
QKLPGNDNTS...	QKLPGNDNSS...


2 Run matrix_generate.py to generate the input matrix from the input file:

python matrix_generate.py 
--aaindex_file aaindex_feature_H1N1/aaindex_feature_H3N2
--seq_file 
--type training/testing/predict
--dir save_dir


3 Run train.py to train a CNN model:

python train.py
--train_data 
--test_data 
--outdir 
--type H1N1/H3N2

4 Run predict.py to predict the relationship between the two sequence 

python predict.py
--test_data 
--seq_file 
--model_path 
--outdir 
--type H1N1/H3N2

Example of the predicted file

seq_1	seq_2	predict
QKLPGNDNST...	QKLPGNDNSS...	0
-KLPGNDNS...	-KLPGNDNT...	1
QKLPGIDNSN...	QKLPGIDNSS...	0
QKLPGNDNTS...	QKLPGNDNSS...	1



Jing Meng<br>

jing.mengrabbit@outlook.com<br>
