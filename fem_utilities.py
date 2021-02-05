import array
import csv
import re
import math
from typing import List, Any

import numpy as np
from copy import deepcopy


def make_nodes(length, breadth, nodes_in_x, nodes_in_y):
    x = np.linspace(0, length, nodes_in_x)
    y = np.linspace(0, breadth, nodes_in_y)
    xx, yy = np.meshgrid(x, y)
    no_of_nodes = len(x) * len(y)
    arr = np.zeros((no_of_nodes + 1, 3))
    count = 1
    for j in range(len(xx)):
        for k in range(len(xx[0])):
            arr[count][0] = int(count)
            arr[count][1] = xx[j][k]
            arr[count][2] = yy[j][k]
            count += 1
    arr = arr[1:].tolist()
    arr = [[int(x[0]), x[1], x[2]] for x in arr]
    return arr


# x = make_nodes(10, 20, 3, 5)
# print(x)


def connectivity(length, breadth, nodes_in_x, nodes_in_y):
    x = np.linspace(0, length, nodes_in_x)
    y = np.linspace(0, breadth, nodes_in_y)
    xx, yy = np.meshgrid(x, y)
    m = len(xx) - 1
    n = len(xx[0]) - 1
    no_of_nodes = len(x) * len(y)
    NElements = m * n
    g = np.reshape(np.arange(1, no_of_nodes + 1), (len(xx), len(xx[0])))
    eles = np.zeros((m * n + 1, 5), dtype=int)
    cnt = 1
    for i in range(len(g) - 1):
        for j in range(len(g[0]) - 1):
            eles[cnt][0] = int(cnt)
            eles[cnt][1] = g[i][j]
            eles[cnt][2] = g[i][j + 1]
            eles[cnt][3] = g[i + 1][j + 1]
            eles[cnt][4] = g[i + 1][j]
            cnt += 1
    eles = eles[1:].tolist()
    eles = [[int(x[0]), int(x[1]), int(x[2]), int(x[3]), int(x[4])] for x in eles]
    return eles, NElements


# c, n= connectivity(10, 20, 3, 5)
# print(c)


def make_BC(Node_Data):
    NDOF = 2 * len(Node_Data)
    BC = [[0 for row in range(1)] for col in range(NDOF)]  # required matrix for BCs
    for x in Node_Data:
        if x[1] == 0:
            BC[2 * (x[0] - 1)][0] = 1
            BC[2 * (x[0] - 1) + 1][0] = 1

    return BC


# y = make_BC(x)
# print(y)
def apply_loads(Node_Data, length):
    NNodes = len(Node_Data)
    NDOF = 2 * NNodes
    l = 0
    for x_ in Node_Data:
        if x_[1] == 0:
            l += 1
    print(l)
    F = 10000 / l
    print(F)
    LOADS = [[0 for row in range(1)] for col in range(NDOF)]  # required matrix for loading conditions
    for x in Node_Data:
        if x[1] == length:
            LOADS[2 * (x[0] - 1)][0] = F
    print(LOADS)
    return LOADS


# loads = apply_loads(x, 10)
# print(len(loads))

def matrixmult(A, B):
    A = np.array(A)
    B = np.array(B)
    C = np.matmul(A, B)
    return C.tolist()


# for adding two matrices
def matrixsum(A, B):
    A = np.array(A)
    B = np.array(B)
    C = A + B
    return C.tolist()


# scalar multiplication with a matrix
def matrixscalarmult(A, b):
    A = np.array(A)
    C = b * A
    return C.tolist()


# for taking a transpose
def transpose(A):
    A = np.array(A)
    A = A.T
    return A.tolist()


def calculate_D(Analysis_Type, E, Poisson):
    if Analysis_Type == 1:
        D = [[E / (1 - Poisson * Poisson), (Poisson * E) / (1 - Poisson * Poisson), 0],
             [(Poisson * E) / (1 - Poisson * Poisson), E / (1 - Poisson * Poisson), 0],
             [0, 0, ((1 - Poisson) / 2 * E) / (1 - Poisson * Poisson)]]
    else:
        D = [
            [E / ((1 + Poisson) * (1 - 2 * Poisson)) * (1 - Poisson), Poisson * E / ((1 + Poisson) * (1 - 2 * Poisson)),
             0],
            [Poisson * (E / ((1 + Poisson) * (1 - 2 * Poisson))),
             E / ((1 + Poisson) * (1 - 2 * Poisson)) * (1 - Poisson), 0],
            [0, 0, E / ((1 + Poisson) * (1 - 2 * Poisson)) * ((1 - 2 * Poisson) / 2)]]

    return D


# print(calculate_D(1, 210e9, 0.3))
def DetJacobian(s, t, X, Y):
    J = [[0, 0], [0, 0]]
    J[0][0] = (t - 1) * X[0][0] + (1 - t) * X[0][1] + (1 + t) * X[0][2] - (1 + t) * X[0][3]
    J[0][1] = (t - 1) * Y[0][0] + (1 - t) * Y[0][1] + (1 + t) * Y[0][2] - (1 + t) * Y[0][3]
    J[1][0] = (s - 1) * X[0][0] - (1 + s) * X[0][1] + (1 + s) * X[0][2] + (1 - s) * X[0][3]
    J[1][1] = (s - 1) * Y[0][0] - (1 + s) * Y[0][1] + (1 + s) * Y[0][2] + (1 - s) * Y[0][3]
    return (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / 16.0


def BMatrix(s, t, X, Y):
    B = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    a = 0.25 * (Y[0][0] * (s - 1) + Y[0][1] * (-1 - s) + Y[0][2] * (1 + s) + Y[0][3] * (1 - s))
    b = 0.25 * (Y[0][0] * (t - 1) + Y[0][1] * (1 - t) + Y[0][2] * (1 + t) + Y[0][3] * (-1 - t))
    c = 0.25 * (X[0][0] * (t - 1) + X[0][1] * (1 - t) + X[0][2] * (1 + t) + X[0][3] * (-1 - t))
    d = 0.25 * (X[0][0] * (s - 1) + X[0][1] * (-1 - s) + X[0][2] * (1 + s) + X[0][3] * (1 - s))

    B[0][0] = (a * 0.25 * (t - 1) - b * 0.25 * (s - 1)) / DetJacobian(s, t, X, Y)
    B[1][0] = 0
    B[2][0] = (c * 0.25 * (s - 1) - d * 0.25 * (t - 1)) / DetJacobian(s, t, X, Y)

    B[0][1] = 0
    B[1][1] = (c * 0.25 * (s - 1) - d * 0.25 * (t - 1)) / DetJacobian(s, t, X, Y)
    B[2][1] = (a * 0.25 * (t - 1) - b * 0.25 * (s - 1)) / DetJacobian(s, t, X, Y)

    B[0][2] = (a * 0.25 * (1 - t) - b * 0.25 * (-1 - s)) / DetJacobian(s, t, X, Y)
    B[1][2] = 0
    B[2][2] = (c * 0.25 * (-1 - s) - d * 0.25 * (1 - t)) / DetJacobian(s, t, X, Y)

    B[0][3] = 0
    B[1][3] = (c * 0.25 * (-1 - s) - d * 0.25 * (1 - t)) / DetJacobian(s, t, X, Y)
    B[2][3] = (a * 0.25 * (1 - t) - b * 0.25 * (-1 - s)) / DetJacobian(s, t, X, Y)

    B[0][4] = (a * 0.25 * (1 + t) - b * 0.25 * (1 + s)) / DetJacobian(s, t, X, Y)
    B[1][4] = 0
    B[2][4] = (c * 0.25 * (1 + s) - d * 0.25 * (1 + t)) / DetJacobian(s, t, X, Y)

    B[0][5] = 0
    B[1][5] = (c * 0.25 * (1 + s) - d * 0.25 * (1 + t)) / DetJacobian(s, t, X, Y)
    B[2][5] = (a * 0.25 * (1 + t) - b * 0.25 * (1 + s)) / DetJacobian(s, t, X, Y)

    B[0][6] = (a * 0.25 * (-1 - t) - b * 0.25 * (1 - s)) / DetJacobian(s, t, X, Y)
    B[1][6] = 0
    B[2][6] = (c * 0.25 * (1 - s) - d * 0.25 * (-1 - t)) / DetJacobian(s, t, X, Y)

    B[0][7] = 0
    B[1][7] = (c * 0.25 * (1 - s) - d * 0.25 * (-1 - t)) / DetJacobian(s, t, X, Y)
    B[2][7] = (a * 0.25 * (-1 - t) - b * 0.25 * (1 - s)) / DetJacobian(s, t, X, Y)

    return B


def make_Ke(ls, lt, X, Y, D, Thickness):
    Ke1 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(ls[0], lt[0], X, Y)), matrixmult(D, BMatrix(ls[0], lt[0], X, Y))),
                         DetJacobian(ls[0], lt[0], X, Y)),
        Thickness)
    Ke2 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(ls[1], lt[1], X, Y)), matrixmult(D, BMatrix(ls[1], lt[1], X, Y))),
                         DetJacobian(ls[1], lt[1], X, Y)),
        Thickness)
    Ke3 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(ls[2], lt[2], X, Y)), matrixmult(D, BMatrix(ls[2], lt[2], X, Y))),
                         DetJacobian(ls[2], lt[2], X, Y)),
        Thickness)
    Ke4 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(ls[3], lt[3], X, Y)), matrixmult(D, BMatrix(ls[3], lt[3], X, Y))),
                         DetJacobian(ls[3], lt[3], X, Y)),
        Thickness)

    return matrixsum(matrixsum(matrixsum(Ke1, Ke2), Ke3), Ke4)


# Mass Matrix
def NMatrix(s, t):
    NMatrix = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]

    NMatrix[0][0] = ((1 - s) * (1 - t) / 4)
    NMatrix[1][0] = 0

    NMatrix[0][1] = 0
    NMatrix[1][1] = ((1 - s) * (1 - t) / 4)

    NMatrix[0][2] = ((1 + s) * (1 - t) / 4)
    NMatrix[1][2] = 0

    NMatrix[0][3] = 0
    NMatrix[1][3] = ((1 + s) * (1 - t) / 4)

    NMatrix[0][4] = ((1 + s) * (1 + t) / 4)
    NMatrix[1][4] = 0

    NMatrix[0][5] = 0
    NMatrix[1][5] = ((1 + s) * (1 + t) / 4)

    NMatrix[0][6] = ((1 - s) * (1 + t) / 4)
    NMatrix[1][6] = 0

    NMatrix[0][7] = 0
    NMatrix[1][7] = ((1 - s) * (1 + t) / 4)

    return NMatrix


def mass_matrix(ls, lt, X, Y, Thickness, density):
    MassMatrix1 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(ls[0], lt[0])), NMatrix(ls[0], lt[0])),
                                          DetJacobian(ls[0], lt[0], X, Y)),
                         Thickness), density)
    MassMatrix2 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(ls[1], lt[1])), NMatrix(ls[1], lt[1])),
                                          DetJacobian(ls[1], lt[1], X, Y)),
                         Thickness), density)
    MassMatrix3 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(ls[2], lt[2])), NMatrix(ls[2], lt[2])),
                                          DetJacobian(ls[2], lt[2], X, Y)),
                         Thickness), density)
    MassMatrix4 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(ls[3], lt[3])), NMatrix(ls[3], lt[3])),
                                          DetJacobian(ls[3], lt[3], X, Y)),
                         Thickness), density)

    MassMatrixelement = matrixsum(matrixsum(matrixsum(MassMatrix1, MassMatrix2), MassMatrix3), MassMatrix4)

    # *************************************************************************************
    # GLOBAL STIFFNESS ANS MASS MATRIX ASSEMBLY
    # *************************************************************************************
    # NDOF = 2 * len(Node_Data)
    # KG = [[0 for row in range(NDOF)] for col in range(NDOF)]
    # Mmatrix = [[0 for row in range(NDOF)] for col in range(NDOF)]
    # Ke = make_Ke(ls, lt, X, Y, Thickness)
    # for r in range(0, 4):
    #     for t in range(0, 2):
    #         for z in range(0, 4):
    #             for m in range(0, 2):
    #                 rw = NDOFN * (EL_Data[i - 1][r + 1] - 1) + t
    #                 cl = NDOFN * (EL_Data[i - 1][z + 1] - 1) + m
    #                 KG[rw][cl] = KG[rw][cl] + Ke[NDOFN * r + t][NDOFN * z + m]
    #                 Mmatrix[rw][cl] = Mmatrix[rw][cl] + MassMatrixelement[NDOFN * r + t][NDOFN * z + m]

    return MassMatrixelement


def remove_force_from_BC(LOADS, BC, NDOF):
    # LOADS_total=np.add(matrixmult(Mmatrix,Acceleration),LOADS)
    LOADS_total = LOADS
    # --------------removing the body forces for the nodes with BCs
    index = 0
    while index < NDOF:
        if BC[index][0] == 1:
            LOADS_total[index][0] = 0
        else:
            LOADS_total[index][0] = LOADS_total[index][0]
        index = index + 1
    # print(LOADS_total)

    return LOADS_total


def apply_BC(BC, KG, NDOF):
    # APPLYING BOUNDARY CONDITIONS
    # *************************************************************************************
    # --------------Filling Global AUXILIAR Stiffness Matrix with zeros
    KG2 = [[0 for row in range(NDOF)] for col in range(NDOF)]
    index = 0
    while index < NDOF:
        if BC[index][0] == 1:
            for r in range(0, NDOF):
                if r == index:
                    KG2[index][r] = 1
                else:
                    KG2[index][r] = 0
                    KG2[r][index] = 0
        else:
            for r in range(index, NDOF):
                KG2[index][r] = KG[index][r]
            for t in range(index, NDOF):
                KG2[t][index] = KG[t][index]
        index = index + 1
    # print (KG2)
    return KG2


def compute_displacement(KG2, LOADS_total):
    KGINV = np.linalg.pinv(KG2)
    displacement = matrixmult(KGINV, LOADS_total)

    return displacement


def max_disp(disp):
    return max(disp[:])[0]


def calc_average(disp, Node_Data, length):
    d = 0
    t = 0
    for x in Node_Data:
        if x[1] == length:
            idx = x[0] - 1
            d += disp[2 * idx][0]
            t += 1
    avg = d / t
    return avg
