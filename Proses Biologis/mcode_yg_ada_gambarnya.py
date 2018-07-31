# Author: True Price <jtprice@cs.unc.edu>

import networkx as nx
import sys
from collections import defaultdict
import numpy as np
import math
from matplotlib.pylab import show, cm, axis
import matplotlib.pyplot as plt

def draw_graph(d, graph, clusters, **kwargs):
    """
    Visualize the clustering

    :param matrix: The unprocessed adjacency matrix
    :param clusters: list of tuples containing clusters as returned
                     by 'get_clusters'
    :param kwargs: Additional keyword arguments to be passed to
                   networkx.draw_networkx
    """
    # make a networkx graph from the adjacency matrix

    # map node to cluster id for colors
    # cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}
    num=0
    color=[]
    for i in graph:
        num=d[i]
        color.append(num)

    colors=np.array(color)
    # if colormap not specified in kwargs, use a default
    if not kwargs.get("cmap", False):
        kwargs["cmap"] = cm.tab20

    # draw
    nx.draw_networkx(graph, node_color=colors, **kwargs)
    ax = plt.gca()
    axis("off")
    plt.show()
    #
    # show(block=False)

WEIGHT_THRESHOLD = 0.2

##
WEIGHT_THRESHOLD = 1 - WEIGHT_THRESHOLD

def algo_mcode(filename):
  edges = defaultdict(set) # node id => neighboring node ids
  dic={}
  G=nx.Graph()

  # read in graph
  with open(filename, 'r') as f:
    for line in f:
      a,b = line.split()[:2]
      edges[a].add(b)
      edges[b].add(a)
      dic[a]=0
      dic[b]=0
      G.add_edge(a,b)
  print >> sys.stderr, 'graph loaded; %i nodes' % (len(edges),)
  
  # Stage 1: Vertex Weighting
  print >> sys.stderr, 'vertex weighting...'
  weights = dict((v,1.) for v in edges)
  for i,v in enumerate(edges):
    if i % 1000 == 0: print >> sys.stderr, i
    neighborhood = set((v,)) | edges[v]
    # if node has only one neighbor, we know everything we need to know
    if len(neighborhood) <= 2: continue

    # see if larger k-cores exist
    k = 3 # highest valid k-core
    while neighborhood:
      k_core = neighborhood.copy()
      invalid_nodes = True
      while invalid_nodes and neighborhood:
        invalid_nodes = set(
          n for n in neighborhood if len(edges[n] & neighborhood) <= k)
        neighborhood -= invalid_nodes
      k += 1 # on exit, k will be one greater than we want
    # vertex weight = k-core number * density of k-core
    weights[v] = (k-1) * (sum(len(edges[n] & k_core) for n in k_core) /
      (2. * len(k_core)**2))

  # Stage 2: Molecular Complex Prediction
  print >> sys.stderr, 'molecular complex prediction'
  unvisited = set(edges)
  num_clusters = 0
  mc = list()
  colorn=0
  for seed in sorted(weights, key=weights.get, reverse=True):
    if seed not in unvisited: continue

    cluster, frontier = set((seed,)), set((seed,))
    colorn+=1
    dic[seed]=colorn
    w = weights[seed] * WEIGHT_THRESHOLD
    while frontier:
      cluster.update(frontier)
      unvisited -= frontier
      frontier = set(
        n for n in set.union(*(edges[n] for n in frontier)) & unvisited
        if weights[n] > w)

    # haircut: only keep 2-core complexes
    invalid_nodes = True
    while invalid_nodes and cluster:
      invalid_nodes = set(n for n in cluster if len(edges[n] & cluster) < 2)
      cluster -= invalid_nodes

    if cluster:
      # fluff never really seems to improve anything...
      #cluster.update(
      # n for n in set.union(*(edges[c] for c in cluster)) & unvisited
      # if densities[n] > FLUFF_THRESHOLD)

      print(' '.join(cluster))
      num_clusters += 1
      print >> sys.stderr, num_clusters, len(cluster), seed

    for i in cluster:
        dic[i]=num_clusters
    mc.append(cluster)
  return dic, G, mc

def mcode(filename):
    d, g, cluster=algo_mcode(filename)
    a=len(g)
    pos = nx.spring_layout(g, k=2.5/math.sqrt(a))
    draw_graph(d, g, cluster
           , pos=pos, node_size=15, with_labels=False, edge_color="silver")

if __name__ == '__main__':
  mcode(sys.argv[1])
 
