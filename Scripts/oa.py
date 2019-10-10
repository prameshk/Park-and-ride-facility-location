'''

Created on Thurs Oct 10 21:09:46 2019

@author: Pramesh Kumar

This programs uses Outer Approximation to solve the following
Mixed Integer Non-Linear Program related to Discrete network design problem

'''

from gurobipy import * # For solving the optimization problem
import cvxpy as cp
import numpy as np
import math


def NLP(y):
    x = cp.Variable(2)
    objective = cp.Minimize(y[0] + y[1] + cp.square(x[0]) +  cp.square(x[1]))
    constraints = [cp.square(x[0]-2) - x[1] <= 0, x[0] - 2*y[0] >= 0, x[0] - x[1] - 3*(1 - y[0]) <= 0, x[0] - (1 - y[0]) >= 0, x[1] - y[1] >= 0, x[0] + x[1] >= 3*y[0], x[0] >= 0, x[0] <= 4, x[1] >= 0, x[1] <= 4]
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



def OuterApprox(eps, y):
    UB = float("inf")
    LB = -float("inf")
    cuts = []
    tol = float("inf")
    while tol > eps:
        c, u = NLP(y)
        UB = min(UB, u)
        cuts.append(c)
        x, y, l = MasterProblem(cuts)
        LB = max(LB, l)
        tol = UB - LB

    return x, y, UB, LB




print(OuterApprox(0.1, [1, 1]))