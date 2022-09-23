# CNNantigentic

CNNantigentic is based on the CNN to predict the antigenic relationship of influenza A/H3N2 virus. CNNantigenic constructs a spatially oriented representation of the HA1 sequence adapted for the convolutional architecture, which can explore the interactions of the amino acids in the context sequence. Moreover, rather than the redundant amino acid embeddings, CNNantigenic takes into account only physicochemical features determining antigenicity of influenza A/H3N2 virus.Together, CNNantigenic can effectively extract the features in the context sequence from local to global views, and investigate the combinatorial contributions of point mutations in the HA protein to the antigenicity.

## Prerequisites

TensorFlow 2.7.0 
sklearn 0.19.2 
pandas 1.3.5 
numpy 1.21.6 
keras 2.7.0 

## Preview

├─model
└─script
    └─input_matrix_generation

The code of the overall CNNantigenic project includes two parts, in which the **model** stores the prediction model saved from 2006 to 2020 (for example, 2020.h5 represents the model of CNNantigenic that predicts the antigen relationship in 2020 year), and **script** stores the specific script, where **input_ matrix_ generation** is used to convert H3N2 correlation data into the input matrix of the model.

## Usage

There are three steps to predict antigen relationship

1. In the folder **input_ matrix_ generation**. Using files <u>csv_to_5fold.py</u> and <u>csv_to_years.py</u>, the total data can be divided into the separate data required for the 5-cross validation and retrospective testing.And then, we can use files <u>matrix_generation_single.py</u> and <u>matrix_generation_double.py</u> to convert separate data to input_matrix. The difference between <u>matrix_generation_single.py</u> and <u>matrix_generation_double.py</u>  in whether data amplification is performed. Finally, we will get the input matrix of numpy type
2. In the folder **script**. The file of <u>main.py</u> used to define model, training model and test model. When we have prepared the data in step 1, we can use *run_years(years1, years2)* and *run_5folds()* to train the 5-cross validation and retrospective testing respectively. The model will be saved for subsequent testing. We can also *run_ test_Years (years)* to predict the antigen of the year's H3N2 virus by the year's model
3. The file of <u>ROC.py</u> is used to test the specific performance of the saved model, including ROC curve and other indicators

## Example

1. Use <u>csv_to_years.py</u> we can convert  <u>exmple_all_data.csv</u> to  <u>exmple_2020_Relation.csv</u>
2. Use <u>matrix_generation_double.py</u> we can convert <u>exmple_2020_Relation.csv</u> to <u>exmple_2020_double.npy</u>
3. Use the function *run_ test_Years (years)* in the file of <u>main.py</u> we can get the prediction about 2020 year
4. Draw ROC curve by ROC.py

Jing Meng<br>

jing.mengrabbit@outlook.com<br>
