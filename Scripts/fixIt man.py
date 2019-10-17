# -*- coding: utf-8 -*-
"""
Created on Sun October 10 21:09:46 2019

@author: Pramesh Kumar
"""
import math
import time
import heapq
import numpy as np
from scipy import optimize

inputLocation = "Data/"



class Zone:
    def __init__(self, _tmpIn):
        self.zoneId = _tmpIn[0]
        self.lat = 0
        self.lon = 0
        self.destList = []


class Node:
    '''
    This class has attributes associated with any node
    '''
    def __init__(self, _tmpIn):
        self.Id = _tmpIn[0]
        self.lat = 0
        self.lon = 0
        self.outLinks = []
        self.inLinks = []
        self.label = float("inf")
        self.pred = ""
        self.inDegree = 0
        self.outDegree = 0
        self.order = 0 # Topological order
        self.wi = 0.0 # Weight of the node in Dial's algorithm
        self.xi = 0.0 # Toal flow crossing through this node in Dial's algorithm


class Link:
    '''
    This class has attributes associated with any link
    '''
    def __init__(self, _tmpIn):
        self.tailNode = _tmpIn[0]
        self.headNode = _tmpIn[1]
        self.capacity = float(_tmpIn[2]) # veh per hour
        self.length = float(_tmpIn[3]) # Length
        self.fft = float(_tmpIn[4]) # Free flow travel time (min)
        self.beta = float(_tmpIn[6])
        self.alpha = float(_tmpIn[5])
        self.speedLimit = float(_tmpIn[7])
        #self.toll = float(_tmpIn[9])
        #self.linkType = float(_tmpIn[10])
        self.flow = 1.1
        self.cost =  float(_tmpIn[4]) #float(_tmpIn[4])*(1 + float(_tmpIn[5])*math.pow((float(_tmpIn[7])/float(_tmpIn[2])), float(_tmpIn[6])))
        self.logLike = 0.0
        self.reasonable = True # This is for Dial's stochastic loading
        self.wij = 0.0 # Weight in the Dial's algorithm
        self.xij = 0.0 # Total flow on the link for Dial's algorithm
        self.type = 'Road'
        self.park = 0
        self.fare = 0


class ParkRide:
    def __init__(self, _tmpIn):
        self.id = _tmpIn
        self.active = 1.0




class Demand:
    def __init__(self, _tmpIn):
        self.fromZone = _tmpIn[0]
        self.toNode = _tmpIn[1]
        self.demand = float(_tmpIn[2])

def readDemand():
    inFile = open(inputLocation+ "demand.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        tripSet[tmpIn[0], tmpIn[1]] = Demand(tmpIn)
        if tmpIn[0] not in zoneSet:
            zoneSet[tmpIn[0]] = Zone([tmpIn[0]])
        if tmpIn[1] not in zoneSet:
            zoneSet[tmpIn[1]] = Zone([tmpIn[1]])
        if tmpIn[1] not in zoneSet[tmpIn[0]].destList:
            zoneSet[tmpIn[0]].destList.append(tmpIn[1])

    inFile.close()
    print(len(tripSet), "OD pairs")
    print(len(zoneSet), "zones")

def readNetwork():
    inFile = open(inputLocation + "network.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkSet[tmpIn[0], tmpIn[1]] = Link(tmpIn)
        if tmpIn[0] not in nodeSet:
            nodeSet[tmpIn[0]] = Node(tmpIn[0])
        if tmpIn[1] not in nodeSet:
            nodeSet[tmpIn[1]] = Node(tmpIn[1])
        if tmpIn[1] not in nodeSet[tmpIn[0]].outLinks:
            nodeSet[tmpIn[0]].outLinks.append(tmpIn[1])
        if tmpIn[0] not in nodeSet[tmpIn[1]].inLinks:
            nodeSet[tmpIn[1]].inLinks.append(tmpIn[0])

    inFile.close()

def readTransitShortestPath():
    inFile = open(inputLocation + "transitShortestPath.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        TSP[tmpIn[0], tmpIn[1]] = float(tmpIn[2])


def augmentNetwork(pnrNodes):
    OD_pairs = list(tripSet.keys())
    for t in OD_pairs:
        nodeSet[t] = Node([t])
        nodeSet[t].inLinks =  list(set(nodeSet[t].inLinks + [t[1]] + pnrNodes))
        for p in list(set([t[1]] + pnrNodes)):
            nodeSet[p].outLinks.append(t)
            if p == t[1]:
                linkSet[p, t] = Link([p, t, 0, 0, 0.0005, 0.15, 4, 0])
                linkSet[p, t].park = parkingCost
                linkSet[p, t].type = 'Virtual'
            else:
                try:
                    linkSet[p, t] = Link([p, t, 0, 0, TSP[p, t[1]], 0.15, 4, 0])
                except KeyError:
                    linkSet[p, t] = Link([p, t, 0, 0,100000, 0.15, 4, 0])
                linkSet[p, t].type = 'Transit'
                linkSet[p, t].fare = transitFare

    for t in OD_pairs:
        tripSet[t[0], t] = Demand([t[0], t, tripSet[t].demand])
        del tripSet[t]
        zoneSet[t[0]].destList = []

    for t in tripSet:
        zoneSet[t[0]].destList.append(t[1])








    print(len(nodeSet), "nodes")
    print(len(linkSet), "links")
    print(len(tripSet), "modified O-D pairs")

###########################################################################################################################


#############################################################################################################################
#############################################################################################################################

def DijkstraHeap(origin):
    '''
    Calcualtes shortest path from an origin to all other destinations.
    The labels and preds are stored in node instances.
    '''
    for n in nodeSet:
        nodeSet[n].label = float("inf")
        nodeSet[n].pred = ""
    nodeSet[origin].label = 0.0
    nodeSet[origin].pred = "NA"
    #SE = [(0, origin)]
    SE = [origin]
    while SE:
        currentNode = SE.pop(0)#heapq.heappop(SE)[1]
        currentLabel = nodeSet[currentNode].label
        for toNode in nodeSet[currentNode].outLinks:
            link = (currentNode, toNode)
            newNode = toNode
            newPred =  currentNode
            existingLabel = nodeSet[newNode].label
            newLabel = currentLabel + linkSet[link].cost
            if newLabel < existingLabel:
                #if (newLabel, newNode) not in SE:
                if newNode not in SE:
                    SE.append(newNode)
                    #heapq.heappush(SE, (newLabel, newNode))
                nodeSet[newNode].label = newLabel
                nodeSet[newNode].pred = newPred


def updateTravelTime(type_eq= "UE"):
    '''
    This method updates the travel time on the links with the current flow
    '''
    for l in linkSet:
        if type_eq == "UE":
            if linkSet[l].type == 'Road':
                linkSet[l].cost = linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta))

            else:
                linkSet[l].cost = (1 / theta) * math.log(linkSet[l].flow + 2)  + linkSet[l].fft
        else:
            if linkSet[l].type == 'Road':
                derivative  = linkSet[l].fft * linkSet[l].alpha * linkSet[l].beta * math.pow((linkSet[l].flow*1.0 /linkSet[l].capacity), linkSet[l].beta)
                linkSet[l].cost = derivative + linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta))
                #linkSet[l].cost = linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta)) +  linkSet[l].fft* linkSet[l].alpha * linkSet[l].beta * math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta)
                #linkSet[l].fft*(1 + (linkSet[l].alpha + linkSet[l].alpha  * linkSet[l].beta) * math.pow((linkSet[l].flow /linkSet[l].capacity), linkSet[l].beta))
            else:
                linkSet[l].cost = (1 / theta) * math.log(2 + linkSet[l].flow)  + linkSet[l].fft + (1 / (linkSet[l].beta * (2 + linkSet[l].flow)))



def calculateCost():
    '''
    This function calculates the total system travel time spent in the network
    :return:
    '''
    cost = 0
    for l in linkSet:
        if linkSet[l].type == 'Road':
            cost = cost + linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta)) * linkSet[l].flow
    return cost

def loadAON():
    '''
    This method produces auxiliary flows for all or nothing loading.
    '''
    x_bar = {l: 0.0 for l in linkSet}
    SPTT = 0.0
    for r in originZones:
        DijkstraHeap(r)
        for s in zoneSet[r].destList:
            try:
                dem = tripSet[r, s].demand
            except KeyError:
                dem = 0.0
            SPTT = SPTT + nodeSet[s].label * dem
            if r != s:
                for spLink in tracePreds(s):
                    x_bar[spLink] = x_bar[spLink] + dem
    return SPTT, x_bar


def tracePreds(dest):
    '''
    This method traverses predecessor nodes in order to create a shortest path
    '''
    prevNode = nodeSet[dest].pred
    spLinks = []
    while nodeSet[dest].pred != "NA":
        spLinks.append((prevNode, dest))
        dest = prevNode
        prevNode = nodeSet[dest].pred
    return spLinks

def findReasonableLinks():
    for l in linkSet:
        if nodeSet[l[1]].label > nodeSet[l[0]].label:
            linkSet[l].reasonable = True
        else:
            linkSet[l].reasonable = False

def computeLogLikelihood():
    '''
    This method computes link likelihood for the Dial's algorithm
    '''
    for l in linkSet:
        if linkSet[l].reasonable == True: # If reasonable link
            linkSet[l].logLike = math.exp(nodeSet[l[1]].label - nodeSet[l[0]].label - linkSet[l].cost)


def topologicalOrdering():
    '''
    * Assigns topological order to the nodes based on the inDegree of the node
    * Note that it only considers reasonable links, otherwise graph will be acyclic
    '''
    for e in linkSet:
        if linkSet[e].reasonable == True:
                nodeSet[e[1]].inDegree = nodeSet[e[1]].inDegree + 1
    order = 0
    SEL = [k for k in nodeSet if nodeSet[k].inDegree == 0]
    while SEL:
        i = SEL.pop(0)
        order = order + 1
        nodeSet[i].order = order
        for j in nodeSet[i].outLinks:
            if linkSet[i, j].reasonable == True:
                nodeSet[j].inDegree = nodeSet[j].inDegree - 1
                if nodeSet[j].inDegree == 0:
                    SEL.append(j)
    if order < len(nodeSet):
        print("the network has cycle(s)")

def resetDialAttributes():
    for n in nodeSet:
        nodeSet[n].inDegree = 0
        nodeSet[n].outDegree = 0
        nodeSet[n].order = 0
        nodeSet[n].wi = 0.0
        nodeSet[n].xi = 0.0
    for l in linkSet:
        linkSet[l].logLike = 0.0
        linkSet[l].reasonable = True
        linkSet[l].wij = 0.0
        linkSet[l].xij = 0.0



def DialLoad():
    '''
    This method runs the Dial's algorithm and prepare a stochastic loading.
    '''
    resetDialAttributes()
    x_bar = {l: 1.1 for l in linkSet}
    for r in originZones:
        DijkstraHeap(r)
        findReasonableLinks()
        topologicalOrdering()
        computeLogLikelihood()



        '''
        Assigning weights to nodes and links
        '''
        order = 1
        while (order <= len(nodeSet)):
            i = [k for k in nodeSet if nodeSet[k].order == order][0] # Node with order no equal to current order
            if order == 1:
                nodeSet[i].wi = 1.0
            else:
                nodeSet[i].wi = sum([linkSet[k, i].wij for k in nodeSet[i].inLinks if linkSet[k, i].reasonable == True])
            for j in nodeSet[i].outLinks:
                if linkSet[i, j].reasonable == True:
                    linkSet[i, j].wij = nodeSet[i].wi*linkSet[i, j].logLike
            order = order + 1
        '''
        Assigning load to nodes and links
        '''
        order = len(nodeSet) # The loading works in reverse direction
        while (order >= 1):
            j = [k for k in nodeSet if nodeSet[k].order == order][0]  # Node with order no equal to current order
            try:
                dem = tripSet[r, j].demand
            except KeyError:
                dem = 0.0
            nodeSet[j].xj = dem + sum([linkSet[j, k].xij for k in nodeSet[j].outLinks if linkSet[j, k].reasonable == True])
            for i in nodeSet[j].inLinks:
                if linkSet[i, j].reasonable == True:
                    linkSet[i, j].xij = nodeSet[j].xj * (linkSet[i, j].wij / nodeSet[j].wi)
            order = order - 1
        for l in linkSet:
            if linkSet[l].reasonable == True:
                x_bar[l] = round(x_bar[l] + linkSet[l].xij, 3)

    return x_bar



def assignment(loading, algorithm,type_eq, accuracy = 0.01, maxIter=100):
    '''
    * Performs traffic assignment
    * Type is either deterministic or stochastic
    * Algorithm can be MSA or FW
    * Accuracy to be given for convergence
    * maxIter to stop if not converged
    '''
    it = 1
    gap = float("inf")
    x_bar = {l: 1.1 for l in linkSet}
    startP = time.time()
    while gap > accuracy:
        if algorithm == "MSA" or it < 2:
            alpha = (1.0/it)
        elif algorithm == "FW":
            alpha = findAlpha(x_bar)
        prevLinkFlow = np.array([linkSet[l].flow for l in linkSet])
        for l in linkSet:
            linkSet[l].flow = alpha*x_bar[l] + (1-alpha)*linkSet[l].flow

        updateTravelTime(type_eq)


        if loading == "deterministic":
            SPTT, x_bar = loadAON()
            TSTT = round(sum([linkSet[a].flow * linkSet[a].cost for a in linkSet]), 3)
            SPTT = round(SPTT, 3)
            if it == 1:
                gap = gap + float("inf")
            else:
                gap = round(abs((TSTT / SPTT) - 1), 5)
            #print(TSTT, SPTT, gap)
        elif loading == "stochastic":
            x_bar = DialLoad()
            currentLinkFlow = np.array([linkSet[l].flow for l in linkSet])
            change = (prevLinkFlow -currentLinkFlow)
            if it < 3:
                gap = gap + float("inf")
            else:
                gap = round(np.linalg.norm(np.divide(change, prevLinkFlow,  out=np.zeros_like(change), where=prevLinkFlow!=0)), 2)

        else:
            print("Terminating the program.....")
            print("The loading ", loading, " is unknown")

        it = it + 1
        if it > maxIter:
            print("The assignment did not converge with the desired gap and max iterations are reached")
            print("current gap ", gap)

    if it < maxIter:
        print("Assignment took", time.time() - startP, " seconds")
        print("assignment converged in ", it, " iterations")

###########################################################################################################################
###########################################################################################################################
tripSet = {}
zoneSet = {}
linkSet = {}
nodeSet = {}
TSP = {}

parkingCost = 15
transitFare = 3
theta = 1.0
readDemand()
readNetwork()
readTransitShortestPath()
augmentNetwork(['4'])
originZones = set([k[0] for k in tripSet])
assignment("deterministic", "MSA", "SO", accuracy=0.01, maxIter=100000)

print(calculateCost())




'''
# Writing the results


for l in linkSet:
    if not isinstance(l[1], tuple):
        print(l)


flow = {l:linkSet[l].flow for l in linkSet if not isinstance(l[1], tuple)}

outFile = open("pnrflow.dat", "w")  # IVT, WT, WK, TR
tmpOut = "tailNode\theadNode\tUE_flow"
outFile.write(tmpOut + "\n")
for i in flow:
    tmpOut = i[0] + "\t" + i[1] + "\t" + str(flow[i])
    outFile.write(tmpOut + "\n")
outFile.close()

'''