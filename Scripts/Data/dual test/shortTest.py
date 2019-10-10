# -*- coding: utf-8 -*-
"""
Created on Sun May 28 21:09:46 2017

@author: Pramesh Kumar
"""
import math
import time
import heapq
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
        self.flow = 0.0
        self.cost =  float(_tmpIn[4]) #float(_tmpIn[4])*(1 + float(_tmpIn[5])*math.pow((float(_tmpIn[7])/float(_tmpIn[2])), float(_tmpIn[6])))
        self.linkLike = 0.0


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
        if tmpIn[1] not in zoneSet[tmpIn[1]].destList:
            zoneSet[tmpIn[1]].destList.append(tmpIn[0])

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
    print(len(nodeSet), "nodes")
    print(len(linkSet), "links")



###########################################################################################################################

readStart = time.time()

tripSet = {}
zoneSet = {}
linkSet = {}
nodeSet ={}



readDemand()
readNetwork()

originZones = set([k[0] for k in tripSet])
print("Reading the network data took", round(time.time() - readStart, 2), "secs")

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
    SE = [(0, origin)]
    while SE:
        currentNode = heapq.heappop(SE)[1]
        currentLabel = nodeSet[currentNode].label
        for toNode in nodeSet[currentNode].outLinks:
            link = (currentNode, toNode)
            newNode = toNode
            newPred =  currentNode
            existingLabel = nodeSet[newNode].label
            newLabel = currentLabel + linkSet[link].cost
            if newLabel < existingLabel:
                heapq.heappush(SE, (newLabel, newNode))
                nodeSet[newNode].label = newLabel
                nodeSet[newNode].pred = newPred

def updateTravelTime():
    for l in linkSet:
        linkSet[l].cost = linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta))


def findAlpha(x_bar):
    alpha = 0.0
    def df(alpha):
        sum_derivative = 0 ## this line is the derivative of the objective function.
        for l in linkSet:
            sum_derivative = sum_derivative + (x_bar[l] - linkSet[l].flow)*BPR((linkSet[l].flow + alpha*(x_bar[l] - linkSet[l].flow)), linkSet[l].fft, linkSet[l].alpha, linkSet[l].beta, linkSet[l].capacity)
        return sum_derivative
    sol = optimize.root(df, 0)
    return max(0, min(1, sol.x[0]))

def tracePreds(dest):
    prevNode = nodeSet[dest].pred
    spLinks = []
    while nodeSet[dest].pred != "NA":
        spLinks.append((prevNode, dest))
        dest = prevNode
        prevNode = nodeSet[dest].pred
    return spLinks

def loadAON():
    x_bar = {l: 0.0 for l in linkSet}
    SPTT = 0.0
    for r in originZones:
        DijkstraHeap(r)
        for s in zoneSet[r].destList:
            dem = tripSet[r, s].demand
            SPTT = SPTT + nodeSet[s].label * dem
            if r != s:
                for spLink in tracePreds(s):
                    x_bar[spLink] = x_bar[spLink] + dem
    return SPTT, x_bar


def computeLinkLikelihood():
    '''
    This method computes link likelihood for the Dial's algorithm
    '''


def dialLoad():
    for r in originZones:
        DijkstraHeap(r)

    pass


def assignment(loading, algorithm, accuracy = 0.01, maxIter=100):
    '''
    * Performs traffic assignment
    * Type is either deterministic or stochastic
    * Algorithm can be MSA or FW
    * Accuracy to be given for convergence
    * maxIter to stop if not converged
    '''
    it = 1
    gap = float("inf")
    x_bar = {l: 0.0 for l in linkSet}
    startP = time.time()
    while gap > accuracy:
        if algorithm == "MSA" or it < 2:
            alpha = (1.0/it)
        elif algorithm == "FW":
            alpha = findAlpha(x_bar)
        else:
            print("Terminating the program.....")
            print("The solution algorithm ", algorithm, " does not exist!")
        for l in linkSet:
            linkSet[l].flow = alpha*x_bar[l] + (1-alpha)*linkSet[l].flow
        updateTravelTime()
        if loading == "deterministic":
            SPTT, x_bar = loadAON()
        elif loading == "stochastic":
            SPTT, x_bar = loadAON()
        else:
            print("Terminating the program.....")
            print("The loading ", loading, " is unknown")
        TSTT = round(sum([linkSet[a].flow*linkSet[a].cost for a in linkSet]), 3)
        SPTT = round(SPTT, 3)
        gap = round(abs((TSTT / SPTT) - 1), 5)
        #print(TSTT, SPTT, gap)
        if it == 1:
            gap = gap  + float("inf")
        it = it + 1
        if it > maxIter:
            print("The assignment did not converge with the desired gap")
            print("current gap ", gap)
            break
    print("Assignment took", time.time() - startP)
    return "assignment converged in ", it, " iterations"

###########################################################################################################################



#assignment("deterministic", "MSA", accuracy = 0.01, maxIter=100)

linkSet

from gurobipy import *

m = Model()

decVars ={}

constr ={}

for i in linkSet:
    decVars[i] = m.addVar(vtype=GRB.CONTINUOUS, name=str(i),  obj = linkSet[i].fft)

m.update()
constr['1'] = m.addConstr(quicksum([decVars[('1', k)] for k in nodeSet['1'].outLinks]) == 1)
constr['24'] = m.addConstr(quicksum([decVars[(k, '24')] for k in nodeSet['24'].inLinks]) == 1)
constr['100'] = m.addConstr(quicksum([decVars[('24', k)] for k in nodeSet['24'].outLinks]) == 0)
constr['2400'] = m.addConstr(quicksum([decVars[(k, '1')] for k in nodeSet['1'].inLinks]) == 0)
for i in nodeSet:
    if i not in ['1', '24']:
        constr[i] = m.addConstr(quicksum([decVars[(i, k)] for k in nodeSet[i].outLinks]) == quicksum([decVars[(k, i)] for k in nodeSet[i].inLinks]))


m.update()
obj = m.getObjective()
m.setObjective(obj, sense=GRB.MINIMIZE)
m.optimize()

