'''

Created on Thurs Oct 10 21:09:46 2019

@author: Pramesh Kumar

This programs uses Outer Approximation to solve the following
Mixed Integer Non-Linear Program related to Discrete network design problem with SUE

'''

from gurobipy import * # For solving the optimization problem
#import cvxpy as cp
import numpy as np
import math
import time
import heapq
import numpy as np




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
pnrNodes =  ['4', '5', '6', '14', '15', '17', '19', '22'] # ['4', '6', '15', '19']
readDemand()
readNetwork()
readTransitShortestPath()
augmentNetwork(pnrNodes)
originZones = set([k[0] for k in tripSet])
destZones = set([k[1] for k in tripSet])


def calculateCost(flow):
    '''
    This function calculates the total system travel time spent in the network
    '''
    cost = 0
    for l in linkSet:
        if linkSet[l].type == 'Road':
            print(linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((flow[l]*1.0/linkSet[l].capacity), linkSet[l].beta)))
            cost = cost + linkSet[l].fft*(1 + linkSet[l].alpha*math.pow((flow[l]*1.0/linkSet[l].capacity), linkSet[l].beta)) * flow[l]
    return cost



def NLP(y, type = "UE"):
    import PNR_equilibrium as eq
    pnrSol = [k for k in y if y[k] == 1]
    c, u = eq.runAssignment(pnrSol, type)
    eq.resetEverything()
    return c, u


def findFeasibleY(cutY):
    m1 = Model()
    y = {p: m1.addVar(vtype=GRB.BINARY,  name='p', lb = 0, ub = 1) for p in pnrNodes}
    for y_sol in cutY:
        tmp = 0
        print(y_sol)
        for n in y_sol:
            tmp = tmp + (1 - y_sol[n]) * y[n] + y_sol[n] * (1 - y[n])
        m1.addConstr(tmp >= 2)

    m1.addConstr(sum([y[k] for k in pnrNodes]) <= 5)
    m1.update()
    m1.setObjective(0, sense=GRB.MINIMIZE)
    m1.Params.OutputFlag = 0
    m1.optimize()
    try:
        return {k: y[k].x for k in y}
    except AttributeError:
        return {}



def updateTravelTime(flow, type= "UE"):
    '''
    This method updates the travel time on the links with the current flow
    '''
    toReturn = 0
    for l in linkSet:
        if linkSet[l].type == 'Road':
            toReturn = toReturn + linkSet[l].fft * (1 + linkSet[l].alpha*math.pow((flow[l]*1.0/linkSet[l].capacity), linkSet[l].beta))
            if type == "SO":
                toReturn = toReturn +  linkSet[l].fft * (linkSet[l].alpha*linkSet[l].beta *math.pow((flow[l]*1.0/linkSet[l].capacity), linkSet[l].beta))
                #linkSet[l].cost = linkSet[l].cost + linkSet[l].fft * (linkSet[l].alpha* (linkSet[l].beta + 1) * math.pow((linkSet[l].flow*1.0/linkSet[l].capacity), linkSet[l].beta))
        '''
        else:
            if flow[l] == 0:
                toReturn = toReturn + (1 / theta) * math.log(2) + linkSet[l].fft
            else:
                toReturn = toReturn + (1 / theta) * math.log(flow[l])  + linkSet[l].fft
            if type == "SO":
                toReturn = toReturn + (1/theta)
                #linkSet[l].cost = linkSet[l].cost + (1 / theta) * (1 + math.log (linkSet[l].flow)) + linkSet[l].fft
        '''

    return toReturn




def MasterProblem(cutX, cutY, upper, lower):
    m = Model()
    x = {l: m.addVar(vtype=GRB.CONTINUOUS, name=str(l), lb=0) for l in linkSet}
    x_od = {(l, o): m.addVar(vtype=GRB.CONTINUOUS, name=str(l) + str(o), lb=0) for l in linkSet for o in originZones}
    y = {p: m.addVar(vtype=GRB.BINARY,  name=str(p), lb = 0, ub = 1) for p in pnrNodes}
    mu = m.addVar(vtype=GRB.CONTINUOUS, name='mu')
    m.update()

    for c in cutX:
        tmp1 = 0
        tmp2 = 0
        for l in linkSet:
            try:
                flow = c[l]
                if linkSet[l].type == 'Road':
                    fft, cap, alpha, beta = linkSet[l].fft, linkSet[l].capacity,  linkSet[l].alpha, linkSet[l].beta
                    tmp1 = tmp1 + fft * flow * (1 +alpha * np.power(flow / cap, beta)) + fft * (1 + (beta + 1) * alpha  * np.power (flow / cap, beta)) * (x[l] - flow)

                else:
                    tmp2  = tmp2 + flow * ((1 / theta) * math.log(flow) + linkSet[l].fft) + ((1 / theta) * (1 + math.log(flow)) + linkSet[l].fft) * (x[l] - flow)
            except KeyError:
                continue

        m.addConstr(mu >= tmp1 + tmp2)

    # Adding "no-good cuts"

    for y_sol in cutY:
        tmp = 0
        for n in y_sol:
            tmp = tmp + (1 - y_sol[n]) * y[n] + y_sol[n] * (1 - y[n])
        m.addConstr(tmp >= 2)


    m.addConstr(sum([y[k] for k in pnrNodes]) <= 5)

    for o in originZones:
        for n in nodeSet:
            tmp = sum([x_od[(n, k), o] for k in nodeSet[n].outLinks]) - sum([x_od[(k, n), o] for k in nodeSet[n].inLinks])
            #tmp = sum([x_od[k] for k in x_od if k[0][0] == n and k[1] == o]) - sum([x_od[k] for k in x_od if k[0][1] == n and k[1] == o])
            if o == n:
                m.addConstr(tmp == sum([tripSet[t].demand for t in tripSet if t[0] == o]))
            elif (o, n) in tripSet:
                m.addConstr(tmp ==  - tripSet[o, n].demand)
            else:
                m.addConstr(tmp == 0)




    for l in linkSet:
        m.addConstr(x[l] == sum([x_od[k] for k in x_od if k[0] == l]))

    for p in pnrNodes:
        for l in [k for k in linkSet if k[0] == p and linkSet[k].type == 'Transit']:
            m.addConstr(x[l] <= 1000000*y[p])


    m.addConstr(mu >= lower)
    m.addConstr(mu <= upper)


    print("Model defined.. starting to solve that now ...")
    obj = mu
    m.update()
    m.setObjective(obj, sense=GRB.MINIMIZE)
    m.Params.OutputFlag = 0
    m.optimize()
    x_val = {l:round(x[l].x,2)  for l in x}
    y_val = {k:y[k].x for k in y}
    ob_val = tmp1.getValue()
    '''
    m.addConstr(y['4'] == 1)
    m.addConstr(y['6'] == 1)
    m.addConstr(y['15'] == 1)
    m.addConstr(y['19'] == 1)
    '''

    #newObj = updateTravelTime(x_val, type="SO")
    #newObj = sum([x_val[l] * linkSet[l].cost for l in x_val if linkSet[l].cost != float("inf")])
    #p = {l:round(x_od[l].x,2)  for l in x_od}
    m.reset(0)
    return x_val, y_val, ob_val



#c, u = NLP({'4':0, '6':0}, "SO")
#cuts_x = [c]
#cuts_y = [{'4': 0, '6': 0}]
#x, y, l = MasterProblem(cuts_x, cuts_y, 7575042, 4575042)

'''
def MasterProblem(cutX, cutY):
    m = Model()
    x = {}
    x_o = {}
    y = {}

    for l in linkSet:
        x[l] = m.addVar(vtype=GRB.CONTINUOUS, name=str(l), lb=0)

    for o in originZones:
        for l in linkSet:
            x_o[l, o] = m.addVar(vtype=GRB.CONTINUOUS, name=str(l) + str(o), lb=0)

    for p in pnrNodes:
        y[p] =  m.addVar(vtype=GRB.BINARY,  name='p', lb = 0, ub = 1)
        m.update()
    alpha = m.addVar(vtype=GRB.CONTINUOUS, name='alpha')
    m.update()

    for c in cutX:
        tmp = 0
        for l in linkSet:
            try:
                flow = c[l]
            except KeyError:
                flow = 2.0
            if linkSet[l].type == 'Road':
                fft, cap = linkSet[l].fft, c[l]
                tmp = tmp + fft * flow * (1 + 0.15 * np.power(flow / cap, 4)) + fft * (1 + 5 * 0.15 * np.power (flow / cap, 4)) * (x[l] - flow)
            else:
                tmp  = tmp + flow * ((1 / theta) * math.log(flow) + linkSet[l].fft) + ((1 / theta) * (1 + math.log(flow)) + linkSet[l].fft) * (x[l] - flow)
        m.addConstr(alpha >= tmp)

    # Adding "no-good cuts"


    for y_sol in cutY:
        tmp = 0
        for n in y_sol:
            tmp = tmp + (1 - y_sol[n]) * y[n] + y_sol[n] * (1 - y[n])
        m.addConstr(tmp >= 2)

    for o in originZones:
        m.addConstr(sum([x_o[l, o] for l in linkSet if l[0]== o]) == sum([tripSet[t].demand for t in tripSet if t[0] == o]))

    for t in tripSet:
        m.addConstr(sum([x_o[k] for k in x_o if k[0][1] == t[1] and k[1] == t[0]]) == tripSet[t].demand)


    for l in linkSet:
        m.addConstr(x[l] == sum([x_o[k] for k in x_o if k[0] == l]))

    for p in pnrNodes:
        for l in [k for k in linkSet if k[0] == p and linkSet[k].type == 'Transit']:
            m.addConstr(x[l] <= 1000000*y[p])
    print(y)
    print("Model defined.. starting to solve that now ...")
    obj = alpha
    m.update()
    m.setObjective(obj, sense=GRB.MINIMIZE)
    m.Params.OutputFlag = 0
    m.optimize()
    x_val = {l:round(x[l].x,2)  for l in x}
    y_val = {k:y[k].x for k in y}
    ob_val = m.objVal
    m.reset(0)
    return x_val, y_val, ob_val

'''


cut = [{l:2 for l in linkSet}]



#############################################################################################################################
#############################################################################################################################

def OuterApprox(eps, cuts_y, maxIter):
    UB = float("inf")
    LB = -float("inf")
    cuts_x = []
    tol = float("inf")
    y_initial = findFeasibleY(cuts_y)
    if len(y_initial) == 0:
        return 0, 0, 0, 0
    else:
        iter = 0
        while tol > eps:
            c, u = NLP(y_initial, "SO")
            print("current UB", u)
            if u <= UB:
                yopt = y_initial
                xopt = c
                UB = u
            # UB = min(UB, u)
            cuts_x.append(c)
            prevLB = LB
            x, y_initial, l = MasterProblem(cuts_x, cuts_y, 100000000, LB) # y_initial is not initial solution but rather iterative soltion of y
            #cuts_x.append(x)
            LB = max(LB, l)
            tol = abs(UB - LB) / UB
            #tol = (LB - prevLB)/prevLB
            print("inner", UB, LB, tol, x, y_initial)
            if iter > maxIter or prevLB - LB == 0:
                print("The required accuracy could not be reached!")
                break
        return x, y_initial, UB, LB


y_sol = [] #[{p:0 for p in pnrNodes}]

#xopt, yopt, UBopt, LBopt = (OuterApprox(0.01, y_sol, 50))
#c, u = NLP(yopt, "SO")








def ParkAndRideLocationOptimization(cuts_y, eps):
    UB = float("inf")
    LB = -float("inf")
    #cuts_y = [y_sol]
    tol = float("inf")
    while tol > eps:
        x, y, u, l = OuterApprox(0.0110, cuts_y, 50)
        if y == 0:
            return x, y, y_opt
        else:
            c, l = NLP(y, "SO")
            if l >= LB:
                y_opt = y
                LB = l
                cuts_y.append(y)
            else:
                cuts_y.append(y)
            #LB = max(l, LB)
            x, u = NLP(y, "UE")
            UB = min(UB, u)
            tol = UB - LB
            if tol <= 0:
                break
            print("Outer", tol, UB, LB, x, y_opt)
        return x, y, y_opt


x, y, y_opt = ParkAndRideLocationOptimization(y_sol, 0.05)


'''

def NLP(y):
    x = cp.Variable(2)
    objective = cp.Minimize(y[0] + y[1] + cp.square(x[0]) +  cp.square(x[1]))
    constraints = [cp.square(x[0]-2) - x[1] <= 0, x[0] - 2*y[0] >= 0, x[0] - x[1] - 3*(1 - y[0]) <= 0, x[0] - (1 - y[0]) >= 0, x[1] - y[1] >= 0, sum([x[k] for k in range(2)]) >= 3*y[0], x[0] >= 0, x[0] <= 4, x[1] >= 0, x[1] <= 4]
    prob = cp.Problem(objective, constraints)
    result = prob.solve()
    return np.round(x.value, 2), np.round(objective.value, 2)


def MasterProblem(cuts):
    m = Model()
    x = {}
    y = {}
    x[0] = m.addVar(vtype=GRB.CONTINUOUS, name='x0', lb = 0, ub = 4)
    x[1] = m.addVar(vtype=GRB.CONTINUOUS, name='x1', lb = 0, ub = 4)
    y[0] = m.addVar(vtype=GRB.BINARY,  name='y0', lb = 0, ub = 1)
    y[1] = m.addVar(vtype=GRB.BINARY, name='y1', lb = 0, ub = 1)
    alpha = m.addVar(vtype=GRB.CONTINUOUS, name='alpha')
    m.update()
    for c in cuts:
        m.addConstr(alpha >= y[0]  + y[1] + np.square(c[0]) + 2*c[0]*(x[0] - c[0]) + np.square(c[1]) + 2*c[1]*(x[1] - c[1]))
        m.addConstr(np.square(c[0] -2) + 2*(c[0] - 2)*(x[0] - c[0]) - x[1] <= 0)
    m.addConstr(x[0] - 2*y[0] >= 0)
    m.addConstr(x[0] - x[1] - 3*(1 - y[0]) <= 0)
    m.addConstr(x[0] - (1 - y[0]) >= 0)
    m.addConstr(x[1] - y[1] >= 0)
    m.addConstr(x[0] + x[1] >= 3*y[0])
    m.addConstr(y[0] + y[1] >= 1)

    obj = alpha
    m.update()
    #m.write('test.lp')
    m.setObjective(obj, sense=GRB.MINIMIZE)
    m.Params.OutputFlag = 0
    m.optimize()
    toReturn = (np.round([x[k].x for k in [0, 1]], 2), np.round([y[k].x for k in [0, 1]]), m.objVal)
    m.reset(0)
    return toReturn



'''
