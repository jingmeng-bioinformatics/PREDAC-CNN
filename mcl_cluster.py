import pandas as pd
import markov_clustering as mc
import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse



parser = argparse.ArgumentParser()
#parser.add_argument('--seq_dir', required=True)
parser.add_argument('--file', required=True, metavar='file')
parser.add_argument('--dir', required=True, metavar='file')
parser.add_argument('--type', required=True, metavar='file')
args = parser.parse_args()

if args.type=='H1N1':
    inflation =1.43
elif args.type=='H3N2':
    inflation =1.40


a=pd.read_table(args.file,sep=',')

G = nx.from_pandas_edgelist(a,'new_name_1','new_name_2','predict')
a_node=pd.DataFrame(list(G),columns=['name'])
matrix = nx.to_scipy_sparse_matrix(G)
result = mc.run_mcl(matrix, inflation=inflation)
clusters = mc.get_clusters(result)

for j in range(len(clusters)):
    for k in clusters[j]:
        a_node.loc[k,'inflation_'+str(inflation)]=j

a_node.to_csv(args.dir+'/'+args.file[:-4]+'_'+args.type+'_'+str(inflation)+'_MCL_clusters.txt',sep='\t')
