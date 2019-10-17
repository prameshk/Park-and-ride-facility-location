import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx


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
        self.flow = 0.0
        self.pnrflow = 0.0

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

def readFlows():
    inFile = open("flow.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkSet[tmpIn[0], tmpIn[1]].flow = float(tmpIn[2])
    inFile.close()

    inFile = open("pnrflow.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkSet[tmpIn[0], tmpIn[1]].pnrflow = float(tmpIn[2])
    inFile.close()





nodeSet = {}
linkSet = {}
readNodes()
readLinks()
readFlows()






pnrNodes =  ('4', '5', '6', '14', '15', '17', '19', '22')
noPNRnodes = tuple(set([k for k in nodeSet]) - set(pnrNodes))
plt.rcParams['figure.figsize'] = 10, 14
firstLinkSet = []
for k in linkSet:
    if k not in firstLinkSet and (k[1], k[0]) not in firstLinkSet:
        firstLinkSet.append(k)

firstLinkSet = tuple(firstLinkSet)

secondLinkSet = tuple(set([k for k in linkSet]) - set(firstLinkSet))
nodeTypes = ['Nodes with park-and-ride', 'Nodes with no park-and-ride']
nodeColors = ['g', 'b']

# assign each node a type and color via a dictionaries
nodeTypeDict = dict(zip(nodeTypes, [pnrNodes, noPNRnodes]))
nodeColorDict = dict(zip(nodeTypes, nodeColors))
nodePos  = dict(zip(nodeSet,[(nodeSet[n].X,nodeSet[n].Y)
                                        for n in nodeSet]))


scale = max(max([linkSet[l].pnrflow for l in linkSet]), max([linkSet[l].flow for l in firstLinkSet]))
weights = tuple([linkSet[l].flow for l in firstLinkSet])




# generate the graph
g = nx.DiGraph()
g.add_nodes_from(nodeSet)
g.add_edges_from(firstLinkSet)
# create image canvas and axes
fig, ax = plt.subplots(1, figsize=(6,7))
plt.axis('off')
# iterate each nodetype, changing colors and labels of the nodes
for nt in nodeTypes:
    # choose nodes and color for each iteration
    nlist = nodeTypeDict[nt]
    ncolor = nodeColorDict[nt]
    # draw the graph
    nx.draw_networkx_nodes(g,
                           pos=nodePos,
                           nodelist=nlist,
                           ax=ax,
                           node_color=ncolor,
                           label=nt)  # the label for each iteration is
                                      # the node type



ed = nx.draw_networkx_edges(g, nodePos, edgelist=firstLinkSet, edge_color=weights, arrows=False, width = 2, edge_cmap=plt.cm.YlOrRd)
#nx.draw_networkx_edges(g, nodePos, edgelist=secondLinkSet, edge_labels=secondLinkSet, arrows=True, width = 2, alpha = 1)
plt.colorbar(ed, ax=ax)
# here is the problem.  The legend does not inherit the colors.
ax.legend(scatterpoints=1, loc='upper center', bbox_to_anchor=(0.8, 0.6, 0.5, 0.5))
plt.show()









''''
G = nx.DiGraph()

G.add_edges_from(firstLinkSet)
f = plt.figure(1)
ax = f.add_subplot(1,1,1)
plt.axis('off')

pos = {n:(nodeSet[n].X, nodeSet[n].Y) for n in nodeSet}


mcl = nx.draw(G, pos, node_color='b', edge_color=weights, width=5.0, edge_color=weights, edge_cmap=plt.cm.Blues,ax=ax)

plt.legend(loc='center')

f.tight_layout()
plt.show()






G = nx.DiGraph()

G.add_edges_from(firstLinkSet)
plt.axis('off')
# Need to create a layout when doing
# separate calls to draw nodes and edges
pos = {n:(nodeSet[n].X, nodeSet[n].Y) for n in nodeSet}
edge_auto = tuple(l for l in linkSet if linkSet[l].type == 'Auto')
scale = max(max([linkSet[l].pnrflow for l in linkSet]), max([linkSet[l].flow for l in linkSet]))
weights =tuple(linkSet[l].flow/scale for l in linkSet)
print({l:linkSet[l].flow - linkSet[l].pnrflow for l in linkSet})
edge_transit = {l:l for l in linkSet if linkSet[l].type == 'Transit'}

#nx.draw_networkx_nodes(G, pos, node_color = 'lightgrey' , node_size = 500)
#nx.draw_networkx_labels(G, pos)
#nx.draw_networkx_edges(G, pos, edgelist=edge_auto, edge_color=weights, arrows=True, edge_cmap=plt.cm.Blues)
from networkx.drawing.nx_agraph import write_dot
#nx.write_dot(G,'multi.dot')
mcl = nx.draw(G, pos, node_color='b', edge_color=weights, width=10.0, edge_cmap=plt.cm.Blues)

#plt.colorbar(mcl)
#nx.draw_networkx_edges(G, pos, edgelist=edge_transit, edge_color='red', arrows=True)
plt.savefig('plt_flow.png', dpi=None, facecolor='w',
    orientation='portrait', papertype=None, format=None,
    transparent=False, bbox_inches=None, pad_inches=0.1,
    frameon=None)
plt.show()

'''