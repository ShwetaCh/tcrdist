## Prep
# conda create -name tcrdist3
# conda activate tcrdist3 
# pip install tcrdist3==0.2.2
# pip install python-louvain 
# pip install cairosvg

## Way to run
# python3 tcrdist3.py
  
## subpackages
#palmotif for CDR3 logo generation
#tcrsampler for generating V-J gene-matched background receptor sets (non productive from human/mouse)
#pwseqdist for efficient and parallelizable pairwise distance computation

## Run ONE time only to get preformatted back- ground sequence files
#from tcrsampler.setup_db import install_all_next_gen
#install_all_next_gen(dry_run = False)

# Check that file for example 1 is available. If not, download it directly from GitHub  
import os
f = 'clonotypes_minervina.tsv'

#if not os.path.isfile(f):
#       os.system('wget https://raw.githubusercontent.com/kmayerb/tcrdist3_book_chapter/main/data/clonotypes_minervina.tsv')
# Make a folder 'data' where some outputs will be written
path = 'data'
if not os.path.isdir(path):
       os.mkdir(path)

import os
import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep

# Load Pandas DataFrame from SI Table 6 (Minervina and Pogorelyy et al., 2021)
f = 'clonotypes_minervina.tsv'
df = pd.read_csv(f,sep ="\t")

# Subset to top 6 epitopes.
list_of_epitopes = ['A01_TTD', 'A02_YLQ', 'A01_LTD',
                    'B15_NQK', 'A01_FTS', 'A24_NYN']

df = df[ df['epitope'].isin(list_of_epitopes)].\
  reset_index(drop = True)

# Rename columns.
df = df.rename(columns = {
    'cdr3b':'cdr3_b_aa',
    'vb'   :'v_b_gene',
    'jb'   :'j_b_gene',
    'cdr3a':'cdr3_a_aa',
    'va'   :'v_a_gene',
    'ja'   :'j_a_gene',
    'donor':'subject'} )

# Preview the first 2 lines of DataFrame.
print(df.head(2))

# Add *01 allele level designation.
df['v_a_gene'] = df['v_a_gene'].apply(lambda x: f"{x}*01")
df['v_b_gene'] = df['v_b_gene'].apply(lambda x: f"{x}*01")
df['j_a_gene'] = df['j_a_gene'].apply(lambda x: f"{x}*01")
df['j_b_gene'] = df['j_b_gene'].apply(lambda x: f"{x}*01")
df['count'] = 1

## Section 4 initialize a TCRrep instance/ compute pairwise distance metrics
from tcrdist.repertoire import TCRrep
tr = TCRrep(cell_df = df[['subject','epitope','cdr3_a_aa','v_a_gene',
                          'j_a_gene','cdr3_b_aa','v_b_gene','j_b_gene',
                          'category','count','cdr3a_nt','cdr3b_nt']],
            organism = 'human',
            chains = ['alpha','beta'],
            deduplicate = True,
            compute_distances = True)

## to see how the four individual CDR matrices are combined to arrive at a weighted multi-CDR distance.

import numpy as np
dst = np.all(tr.pw_beta == (tr.weights_b['cdr1_b_aa']*tr.pw_cdr1_b_aa +
                      tr.weights_b['cdr2_b_aa']*tr.pw_cdr2_b_aa+
                      tr.weights_b['pmhc_b_aa']*tr.pw_pmhc_b_aa +
                      tr.weights_b['cdr3_b_aa']*tr.pw_cdr3_b_aa))

print(dst)

## Section 5 Distance metrics to networks

from tcrdist.public import _neighbors_fixed_radius
# <edge_threshold> is used to define maximum distance to for a network edge.
edge_threshold = 120

# <tr.pw_alpha_beta> is paired chain TCRdist.
tr.pw_alpha_beta = tr.pw_beta + tr.pw_alpha

# <network> initialize a list to populate with edges between TCRs.
network = list()

for i,n in enumerate(_neighbors_fixed_radius(tr.pw_alpha_beta, edge_threshold)):
        for j in n:
                if i != j:
                        network.append((
                                i,                              # 'node_1' - row index
                                j,                              # 'node_2' - col index
                                (tr.pw_alpha_beta )[i,j],       # 'dist'- gets the distance between TCR(i,j)
                                tr.clone_df['v_b_gene'].iloc[i],
                                tr.clone_df['v_b_gene'].iloc[j],
                                tr.clone_df['cdr3_b_aa'].iloc[i],
                                tr.clone_df['cdr3_b_aa'].iloc[j],
                                tr.clone_df['subject'].iloc[i],
                                tr.clone_df['subject'].iloc[j],
                                tr.clone_df['epitope'].iloc[i],
                                tr.clone_df['epitope'].iloc[j],
                                len(n)-1))                      # 'K_neighbors' - number of neighbors
            
            cols = ['node_1', 'node_2', 'dist', 'v_b_gene_1', 'v_b_gene_2',
        'cdr3_b_aa_1','cdr3_b_aa_2', 'subject_1','subject_2',
        'epitope_1','epitope_2', 'K_neighbors']

# Store the <network> edge list as a DataFrame.
df_net = pd.DataFrame(network, columns = cols)

# Option to write the edge list to a file for use in Gephi, Cytoscape, R, etc.
outfile = os.path.join(path,f"{f}_paired_TCRdist_{edge_threshold}_network.csv")
df_net.to_csv(outfile, sep = ",", index = False)

## Add columns to designate edges that unite nodes that come from different donors (i.e., public).
## We can further identify those edges that link TCRs recognizing the same epitope (i.e., consistent edges)

df_net['public'] = df_net['public'] = df_net.apply(lambda x : x['subject_1'] != x['subject_2'], axis = 1)
df_net['consistent'] = df_net.apply(lambda x : x['epitope_1'] == x['epitope_2'], axis = 1)
df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold

##Section 6 Visualization of TCRdist Networks
import matplotlib.pyplot as plt
import networkx as nx
from tcrdist.html_colors import get_html_colors

# Optionally, one can limit network edges to those formed only between TCRs found in two distinct individuals.
df_net = df_net.query('public == True')

# <G> Initialize a networkx Graph instance from the columns of df_net.
G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                          'target' : df_net['node_2'],
                                          'weight' : df_net['weight']}))

# Assign each node a color based on its epitope annotation.
epitopes = ['A01_TTD', 'A02_YLQ', 'A01_LTD', 'B15_NQK', 'A01_FTS', 'A24_NYN']

# Get the same nunber of colors as unique epitopes
colors = get_html_colors(len(epitopes))

# Construct a dictionary to lookup color by epitope.
color_by_epitope = {epitope: color for epitope,color in zip(epitopes,colors)}

# Assign colors to each node based on its epitope annotation.
node_colors = {node: color_by_epitope.get(epitope) for node, epitope in
                zip(df_net['node_1'],df_net['epitope_1'])}

# Positions for all nodes according to a spring layout.
pos = nx.spring_layout(G, seed=2, k = .15)

# Define aesthetic options
options = {"edgecolors": "tab:gray", "node_size": 30, "alpha": 0.5}

#plt.figure()
nx.draw(G,
        nodelist = G.nodes,
        pos = pos,
        node_color= [node_colors[node] for node in G.nodes],
        **options)

plt.savefig('Net1.pdf')

  
