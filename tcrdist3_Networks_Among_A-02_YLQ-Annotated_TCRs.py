import os
import pandas as pd
import networkx as nx
import community.community_louvain as community_louvain
import matplotlib.pyplot as plt
from tcrdist.public import _neighbors_fixed_radius
from tcrdist.repertoire import TCRrep

path = 'data'
f = 'clonotypes_minervina.tsv'

# Only the A*02 YLQ epitope will be considered
epitopes = ["A02_YLQ"]
edge_threshold = 120

df = pd.read_csv(f, sep = "\t")
df = df[ df['epitope'].\
                isin(epitopes)].\
                reset_index(drop = True)

# Rename columns.
df = df.rename(columns = {
        'cdr3b':'cdr3_b_aa',
        'vb':'v_b_gene',
        'jb': 'j_b_gene',
        'cdr3a':'cdr3_a_aa',
        'va': 'v_a_gene',
        'ja' :'j_a_gene',
        'donor':'subject'} )

# Add *01 allele level designation.
df['v_a_gene'] = df['v_a_gene'].apply(lambda x: f"{x}*01")
df['v_b_gene'] = df['v_b_gene'].apply(lambda x: f"{x}*01")
df['j_a_gene'] = df['j_a_gene'].apply(lambda x: f"{x}*01")
df['b_b_gene'] = df['j_b_gene'].apply(lambda x: f"{x}*01")
df['count'] = 1

# <tr> Initialize TCRrep instance.
tr = TCRrep(cell_df = df[['subject','epitope','cdr3_a_aa',
                          'v_a_gene','j_a_gene','cdr3_b_aa',
                          'v_b_gene','j_b_gene','category',
                          'count','cdr3a_nt','cdr3b_nt']] ,
        organism = 'human',
        chains = ['alpha','beta'],
        deduplicate = True,
        compute_distances = True)

# Identify edges among A*02 YLQ annotated TCRs.
network = list()
for i,n in enumerate(_neighbors_fixed_radius(tr.pw_beta+tr.pw_alpha, edge_threshold)):
        for j in n:
                if i != j:
                        network.append((
                        i,
                        j,
                        (tr.pw_beta + tr.pw_alpha)[i,j]
                ))
cols = ['node_1', 'node_2', 'dist']
df_net = pd.DataFrame(network, columns = cols)
df_net['weight'] = edge_threshold - df_net['dist']

G = nx.from_pandas_edgelist(
   pd.DataFrame({'source' : df_net['node_1'],
                 'target' : df_net['node_2'],
                 'weight' :df_net['weight']}))
partition= community_louvain.best_partition(G)

# Change partition such that cluster Id is in descending order based on community size
partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size,
        range(len(partitions_by_cluster_size)))}
partition = {k:partition_reorder.get(v) for k,v in partition.items()}

from tcrdist.html_colors import get_html_colors
clusters = [i for i in pd.Series(partition.values()).value_counts().index]
colors = get_html_colors(len(clusters))
cluster_to_color = {cluster:color for cluster,color, in zip(clusters,colors)}

options = {"edgecolors": "tab:gray","node_size": 50}
pos = nx.spring_layout(G, seed=2, k = .3)
nx.draw(G,
        nodelist = G.nodes,
        pos = pos,
        node_color=[cluster_to_color.get(partition.get(i)) for i in G.nodes],
        **options)

plt.savefig('Net2.pdf')
