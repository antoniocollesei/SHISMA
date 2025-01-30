from dtw import *
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import importlib
np.int = np.int32
import networkx as nx

data = pd.read_csv("temporal_data_ready_normalized.csv", index_col=0)
# remove rows with NaN values
data = data.dropna(axis=0)

index_names = data.index

filter_thresh = None

if filter_thresh is not None:
  array_2d = data.iloc[:filter_thresh, ].to_numpy()
else:
  array_2d = data.to_numpy()

PPI_graph = pd.read_table('FIsInGene_122921_with_annotations.txt')
G = nx.from_pandas_edgelist(PPI_graph, 'Gene1', 'Gene2', 'Score')

# for loop between all time series, calculate the distance between each pair
# and append the distance to a csv file

# create result dataframe with colnames 'index1', 'index2', 'distance'
distance_df = pd.DataFrame(columns=['index1', 'index2', 'distance'])

gene_names = [x.split('_')[0] for x in index_names]
for i in range(len(array_2d)):
  for j in range(i + 1, len(array_2d)):

    if G.has_edge(gene_names[i], gene_names[j]):
      alignment = dtw(array_2d[i], array_2d[j], keep_internals=True)
      # append the distance to the result dataframe
      distance_df = pd.concat([distance_df, pd.DataFrame([{'index1': index_names[i], 'index2': index_names[j], 'distance': alignment.distance}])], ignore_index=True)
      #distance_df = distance_df.append(, ignore_index=True)
    else:
      # skip the pair if there is no edge between the two genes in the PPI network
      continue
    
    # print the progress of the loop
    if j % 10 == 0:
        print("Progress: ", i, j)

# sort distance_df by distance
distance_df = distance_df.sort_values(by='distance', ascending=True)

distance_df.to_csv('tests/distance_df.csv', index=False)