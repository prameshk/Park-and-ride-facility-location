'''

Created on Thurs Oct 10 21:09:46 2019

@author: Pramesh Kumar

This programs uses Outer Approximation to solve the following
Mixed Integer Non-Linear Program related to Discrete network design problem with SUE

'''

from gurobipy import * # For solving the optimization problem
import cvxpy as cp
import numpy as np
import math
import time
import heapq





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
    This class has attributes associated with any Link
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
                    linkSet[p, t] = Link([p, t, 0, 0, 1000000, 0.15, 4, 0])
                linkSet[p, t].type = 'Transit'
                linkSet[p, t].fare = transitFare

    for t in OD_pairs:
        tripSet[t[0], t] = Demand([t[0], t, tripSet[t].demand])
        del tripSet[t]


    print(len(nodeSet), "nodes")
    print(len(linkSet), "links")
    print(len(tripSet), "modified O-D pairs")





###########################################################################################################################
tripSet = {}
zoneSet = {}
linkSet = {}
nodeSet = {}
TSP = {}

parkingCost = 15
transitFare = 3
theta = 1.0
pnrNodes = ['4', '5', '6', '14', '15', '17', '19', '22']
readDemand()
readNetwork()
#readTransitShortestPath()
#augmentNetwork(pnrNodes)

newTripSet ={}

originZones = set([k[0] for k in tripSet])
destZones = set([k[1] for k in tripSet])





def solveSO():
    x = {l: cp.Variable(1) for l in linkSet}
    x_o = {(l, t): cp.Variable(1) for l in linkSet for t in tripSet}
    tempObj = sum([linkSet[l].fft * (x[l] + 0.15 * cp.power(x[l], 5)  / cp.power (linkSet[l].capacity, 4)) for l in linkSet])
    constraints = []
    for t in tripSet:
        for n in nodeSet:
            tmp = sum([x_o[(n, k), t] for k in nodeSet[n].outLinks]) - sum([x_o[(k, n), t] for k in nodeSet[n].inLinks])
            if n == t[0]:
                constraints.append(tmp == tripSet[t].demand)
            elif n == t[1]:
                constraints.append(tmp == - tripSet[t].demand)
            else:
                constraints.append(tmp == 0)

    for l in linkSet:
        constraints.append(x[l] == sum([x_o[k] for k in x_o if k[0] == l]))

    for k in x_o:
        constraints.append(x_o[k] >= 0)

    for k in x:
        constraints.append(x[k] >= 0)

    objective = cp.Minimize(tempObj)
    prob = cp.Problem(objective, constraints)
    print("total decision vars are ", len(x) + len(x_o), " total constraints are ", len(constraints))
    print("Model defined.. starting to solve that now ...")
    result = prob.solve(solver=cp.SCS)


def solveSO():
    x = {l: cp.Variable(1) for l in linkSet}
    x_o = {(l, o): cp.Variable(1) for l in linkSet for o in originZones}
    tempObj = sum([linkSet[l].fft * (x[l] + 0.15 * (cp.power(x[l], 5)  / cp.power (linkSet[l].capacity, 4))) for l in linkSet])
    constraints = []

    for o in originZones:
        for n in nodeSet:
            tmp = sum([ x_o[(n, k), o] for k in nodeSet[n].outLinks]) - sum([ x_o[(k, n), o] for k in nodeSet[n].inLinks])
            if o == n:
                constraints.append(tmp == sum([tripSet[t].demand for t in tripSet if t[0] == o]))
            elif (o, n) in tripSet:
                constraints.append(tmp == - tripSet[o, n].demand)
            else:
                print(n, o)

    for l in linkSet:
        constraints.append(x[l] == sum([x_o[k] for k in x_o if k[0] == l]))

    for k in x_o:
        constraints.append(x_o[k] >= 0)

    for k in x:
        constraints.append(x[k] >= 0)

    objective = cp.Minimize(tempObj)
    prob = cp.Problem(objective, constraints)
    print("Model defined.. starting to solve that now ...")
    result = prob.solve(solver = cp.SCS)
    print(prob.status)









def solveSUE():
    x = {l: cp.Variable(1) for l in linkSet}
    x_o = {(l, o): cp.Variable(1) for l in linkSet for o in originZones}

    tempObj = 0

    for l in linkSet:
        tempObj = tempObj + linkSet[l].fft * ( x[l] + 0.15* (cp.power(x[l], 5) / cp.power(linkSet[l].capacity, 4))) #cp.prod(xij[l], linkSet[l].fft*(1 + linkSet[l].alpha*np.power(xij[l] / linkSet[l].capacity, linkSet[l].beta)))

    constraints = []




    for o in originZones:
        for n in nodeSet:
            tmp = sum([x_o[k] for k in x_o if k[0][0] == n and k[1] == o]) - sum(
                [x_o[k] for k in x_o if k[0][1] == n and k[1] == o])
            print(tmp)
            if o == n:
                constraints.append(tmp== sum([tripSet[t].demand for t in tripSet if t[0] == o]))
            elif (o, n) in tripSet:
                constraints.append(
                    tmp ==  - tripSet[o, n].demand)
            else:
                constraints.append(tmp == 0)

            '''
            if o == n:
                constraints.append(sum([xij_o[l] for l in xij_o if l[1] == o and l[0][0] == n]) - sum([xij_o[l] for l in xij_o if l[1] == o and l[0][1] == n]) == sum([tripSet[t].demand for t in tripSet if t[0] == o]))
            elif n in destZones:
                constraints.append(sum([xij_o[l] for l in xij_o if l[1] == o and l[0][0] == n]) - sum(
                    [xij_o[l] for l in xij_o if l[1] == o and l[0][1] == n]) == tripSet[o, n].demand)
            else:
                print("o yeah")
            '''

    '''
    for t in tripSet:
        constraints.append(sum([xij_o[k] for k in xij_o if k[0][1] == t[1] and k[1] == t[0]]) == tripSet[t].demand)
    '''

    for l in linkSet:
        constraints.append(x[l] == sum([x_o[k] for k in x_o if k[0] == l]))

    for k in x_o:
        constraints.append(x_o[k] >= 0)
    '''
    for p in pnrNodes:
        for l in [k for k in linkSet if k[0] == p and linkSet[k].type == 'Transit']:
            constraints.append(x[l] <= 1000000*y[p])
    '''


    objective = cp.Minimize(tempObj)
    prob = cp.Problem(objective, constraints)
    print("Model defined.. starting to solve that now ...")
    result = prob.solve()






def solveSUE():
    xij = {l: cp.Variable(1) for l in linkSet}
    xij_o = {(l, o): cp.Variable(1) for l in linkSet for o in originZones}

    tempObj = 0

    for l in linkSet:
        tempObj = tempObj + linkSet[l].fft * xij[l] + 0.15* (cp.power(xij[l], 5) / cp.power(linkSet[l].capacity, 4)) #cp.prod(xij[l], linkSet[l].fft*(1 + linkSet[l].alpha*np.power(xij[l] / linkSet[l].capacity, linkSet[l].beta)))

    constraints = []

    for o in originZones:
        for n in nodeSet:
            #constraints.append(sum([xij_o[l] for l in xij_o if l[1] == o and l[0][0] == n]) - sum([xij_o[l] for l in xij_o if l[1] == o and l[0][1] == n]) == sum([tripSet[t].demand for t in tripSet if t[0] == o]) - sum([tripSet[n, d].demand for d in nodeSet[n].inLinks]))
            if o == n:
                constraints.append(sum([xij_o[l] for l in xij_o if l[1] == o and l[0][0] == n]) - sum([xij_o[l] for l in xij_o if l[1] == o and l[0][1] == n]) == sum([tripSet[t].demand for t in tripSet if t[0] == o]))
            else:
                constraints.append(sum([xij_o[l] for l in xij_o if l[1] == o and l[0][0] == n]) - sum(
                    [xij_o[l] for l in xij_o if l[1] == o and l[0][1] == n]) ==  - tripSet[o, n].demand)

    '''
    for t in tripSet:
        constraints.append(sum([xij_o[k] for k in xij_o if k[0][1] == t[1] and k[1] == t[0]]) == tripSet[t].demand)
    '''

    for l in linkSet:
        constraints.append(xij[l] == sum([xij_o[k] for k in xij_o if k[0] == l]))

    for k in xij_o:
        constraints.append(xij_o[k] >= 0)
    '''
    for p in pnrNodes:
        for l in [k for k in linkSet if k[0] == p and linkSet[k].type == 'Transit']:
            constraints.append(x[l] <= 1000000*y[p])
    '''


    objective = cp.Minimize(tempObj)
    prob = cp.Problem(objective, constraints)
    print("Model defined.. starting to solve that now ...")
    result = prob.solve()



