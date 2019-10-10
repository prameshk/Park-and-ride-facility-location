import networkx as nx
import matplotlib.pyplot as plt


class Node:
    '''
    This class has attributes associated with any node
    '''
    def __init__(self, _tmpIn):
        self.Id = _tmpIn[0]
        self.X = float(_tmpIn[1])
        self.Y = float(_tmpIn[2])


class Link:
    '''
    This class has attributes associated with any link
    '''
    def __init__(self, _tmpIn):
        self.tailNode = _tmpIn[0]
        self.headNode = _tmpIn[1]
        self.type = _tmpIn[2]

def readNodes():
    inFile = open("nodes.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        nodeSet[tmpIn[0]] = Node(tmpIn)
    inFile.close()


def readLinks():
    inFile = open("network.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkSet[tmpIn[0], tmpIn[1]] = Link(tmpIn)

    inFile.close()
    print(len(nodeSet), "nodes")
    print(len(linkSet), "links")

plt.rcParams['figure.figsize'] = 10, 14

nodeSet = {}
linkSet = {}
readNodes()
readLinks()

G = nx.MultiGraph()

G.add_edges_from(linkSet)
plt.axis('off')
# Need to create a layout when doing
# separate calls to draw nodes and edges
pos = {n:(nodeSet[n].X, nodeSet[n].Y) for n in nodeSet}
edge_auto = {l:l for l in linkSet if linkSet[l].type == 'Auto'}
edge_transit = {l:l for l in linkSet if linkSet[l].type == 'Transit'}
nx.draw_networkx_nodes(G, pos, node_color = 'lightgrey' , node_size = 500)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edgelist=edge_auto, edge_color='black', arrows=True)
nx.draw_networkx_edges(G, pos, edgelist=edge_transit, edge_color='red', arrows=True)
plt.savefig('plt.png', dpi=None, facecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches=None, pad_inches=0.1,
    frameon=None)
plt.show()

