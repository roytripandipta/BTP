import array
import csv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# Note all units are in SI
Analysis_Type=1  # 1 - plain stress, 2 - plain strain
E = 210e9  # Young modulus
Poisson = 0.3
density = 7.8e-6
Thickness = 1
AccelX = 0
AccelY = 0
NElements = 8
NNodes = 15
NDOF = 2*NNodes
NDOFN = 2
NBC = 5
NLOAD = 5
length = 10
breadth = 20
# make a 2d matrix for material data
Material_data=[[Analysis_Type,E,Poisson,Thickness]]

pi=np.arccos(-1)
if Analysis_Type == 1:
    print('Type of Analysis: Plane Stress')
else:
    print('Type of Analysis: Plane Strain')

# nodal data
NNodes = 0
with open('/content/Nodes_p.csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        NNodes += 1
    f.close()

with open('/content/Nodes_p.csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    Node_Data = [[0 for row in range(3)] for col in range(NNodes)]
    i = 0
    for row in csv_reader:
        Node_Data[i][0] = int(row[0])
        Node_Data[i][1] = float(row[1])
        Node_Data[i][2] = float(row[2])
        i += 1
    f.close()
# print(Node_Data)
# print(NNodes)

# for connectivity matrix
NElements = 0
with open('/content/Eles_p (1).csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        NElements += 1
    f.close()
EL_Data = [[0 for row in range(5)] for col in range(NElements)]
with open('/content/Eles_p (1).csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    i = 0
    for row in csv_reader:
        EL_Data[i][0] = int(row[0])
        EL_Data[i][1] = int(row[1])
        EL_Data[i][2] = int(row[2])
        EL_Data[i][3] = int(row[3])
        EL_Data[i][4] = int(row[4])
        i += 1
    f.close()
# print(EL_Data)

NDOF = 2 * NNodes
# the node numbers on which boundary conditions are imposed
# bc = [1, 4, 7, 10, 13]
BC = [[0 for row in range(1)] for col in range(NDOF)] # required matrix for BCs
for x in Node_Data:
    if x[1] == 0:
        BC[2*(x[0]-1)][0] = 1
        BC[2*(x[0]-1) + 1][0] = 1
# print(BC)

NDOF = 2 * NNodes
# loads = [3,6,9,12,15]
F = 10000/81
LOADS = [[0 for row in range(1)] for col in range(NDOF)] # required matrix for loading conditions
"""
for i in range(NDOF):
  if i+1 in loads:
    LOADS[2*i][0] = F
"""
for x in Node_Data:
    if x[1] == length:
        LOADS[2*(x[0]-1)][0] = F
# print(LOADS)

import numpy as np
# for multiplying 2 matrices
def matrixmult (A, B):
    A = np.array(A)
    B = np.array(B)
    C = np.matmul(A,B)
    return C

# for adding two matrices
def matrixsum (A, B):
   A = np.array(A)
   B = np.array(B)
   C = A + B
   return C

# scaler multiplication with a matrix
def matrixscalarmult(A,b):
   A = np.array(A)
   C = b * A
   return C

# for taking a transpose
def transpose(A):
   A = np.array(A)
   return A.T


# Integration Points Coordinates (4 integration points)

s1=-0.577350269
t1=-0.577350269
s2=0.577350269
t2=-0.577350269
s3=0.577350269
t3=0.577350269
s4=-0.577350269
t4=0.577350269


# intitialize with all zeros
# Global Stiffness matrix
KG = [[0 for row in range(NDOF)] for col in range(NDOF)]

# load matrix
LOADS_total= [[0 for row in range(1)] for col in range(NDOF)]

# Global mass matrix
Mmatrix = [[0 for row in range(NDOF)] for col in range(NDOF)]

# Filling results matrix for each element with zero
Disp_elem = [[0 for row in range(1)] for col in range(8)]
StressX = [[0 for row in range(1)] for col in range(4)]
StressY = [[0 for row in range(1)] for col in range(4)]
StressXY = [[0 for row in range(1)] for col in range(4)]
StressPrinc1 = [[0 for row in range(1)] for col in range(4)]
StressPrinc2 = [[0 for row in range(1)] for col in range(4)]
StressVM = [[0 for row in range(1)] for col in range(4)]
StressX_EXTRAP=[[0 for row in range(1)] for col in range(4)]
StressY_EXTRAP=[[0 for row in range(1)] for col in range(4)]
StressXY_EXTRAP=[[0 for row in range(1)] for col in range(4)]
StressPrinc1_EXTRAP=[[0 for row in range(1)] for col in range(4)]
StressPrinc2_EXTRAP=[[0 for row in range(1)] for col in range(4)]
StressVM_EXTRAP=[[0 for row in range(1)] for col in range(4)]
Strains_integrationpoint1=[[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint2=[[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint3=[[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint4=[[0 for row in range(1)] for col in range(3)]
# PARAVIEW_RESULTS=[[None for row in range(72)] for col in range(NNodes)]

#------------------- nodal extrapolation matrix -------------------------
EXTRAP=[[1+0.5*3**0.5,-0.5,1-0.5*3**0.5,-0.5],[-0.5,1+0.5*3**0.5,-0.5,1-0.5*3**0.5],[1-0.5*3**0.5,-0.5,1+0.5*3**0.5,-0.5],[-0.5,1-0.5*3**0.5,-0.5,1+0.5*3**0.5]]

# Defining material behavior - plane stress or plain strain
if Analysis_Type == 1:
    D = [[E / (1 - Poisson * Poisson), (Poisson * E) / (1 - Poisson * Poisson), 0],
         [(Poisson * E) / (1 - Poisson * Poisson), E / (1 - Poisson * Poisson), 0],
         [0, 0, ((1 - Poisson) / 2 * E) / (1 - Poisson * Poisson)]]
else:
    D = [
        [E / ((1 + Poisson) * (1 - 2 * Poisson)) * (1 - Poisson), Poisson * E / ((1 + Poisson) * (1 - 2 * Poisson)), 0],
        [Poisson * (E / ((1 + Poisson) * (1 - 2 * Poisson))), E / ((1 + Poisson) * (1 - 2 * Poisson)) * (1 - Poisson),
         0],
        [0, 0, E / ((1 + Poisson) * (1 - 2 * Poisson)) * ((1 - 2 * Poisson) / 2)]]

# ========================================================================================
#
#                               ELEMENT STIFFNESS MATRIX
#
# ========================================================================================
for i in range(1, NElements + 1):

    XCOORD = [
        [Node_Data[EL_Data[i - 1][1] - 1][1], Node_Data[EL_Data[i - 1][2] - 1][1], Node_Data[EL_Data[i - 1][3] - 1][1],
         Node_Data[EL_Data[i - 1][4] - 1][1]]]
    YCOORD = [
        [Node_Data[EL_Data[i - 1][1] - 1][2], Node_Data[EL_Data[i - 1][2] - 1][2], Node_Data[EL_Data[i - 1][3] - 1][2],
         Node_Data[EL_Data[i - 1][4] - 1][2]]]


    # ------------------------------------- Jacobian Determinant -------------------------------------------------

    def DetJacobian(s, t):
        J = [[0, 0], [0, 0]]
        J[0][0] = (t - 1) * XCOORD[0][0] + (1 - t) * XCOORD[0][1] + (1 + t) * XCOORD[0][2] - (1 + t) * XCOORD[0][3]
        J[0][1] = (t - 1) * YCOORD[0][0] + (1 - t) * YCOORD[0][1] + (1 + t) * YCOORD[0][2] - (1 + t) * YCOORD[0][3]
        J[1][0] = (s - 1) * XCOORD[0][0] - (1 + s) * XCOORD[0][1] + (1 + s) * XCOORD[0][2] + (1 - s) * XCOORD[0][3]
        J[1][1] = (s - 1) * YCOORD[0][0] - (1 + s) * YCOORD[0][1] + (1 + s) * YCOORD[0][2] + (1 - s) * YCOORD[0][3]
        return (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / 16.0


    # ------------------------------------- Matrix B -------------------------------------------------

    def BMatrix(s, t):
        B = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
        a = 0.25 * (YCOORD[0][0] * (s - 1) + YCOORD[0][1] * (-1 - s) + YCOORD[0][2] * (1 + s) + YCOORD[0][3] * (1 - s))
        b = 0.25 * (YCOORD[0][0] * (t - 1) + YCOORD[0][1] * (1 - t) + YCOORD[0][2] * (1 + t) + YCOORD[0][3] * (-1 - t))
        c = 0.25 * (XCOORD[0][0] * (t - 1) + XCOORD[0][1] * (1 - t) + XCOORD[0][2] * (1 + t) + XCOORD[0][3] * (-1 - t))
        d = 0.25 * (XCOORD[0][0] * (s - 1) + XCOORD[0][1] * (-1 - s) + XCOORD[0][2] * (1 + s) + XCOORD[0][3] * (1 - s))

        B[0][0] = ((a * 0.25 * (t - 1) - b * 0.25 * (s - 1))) / DetJacobian(s, t)
        B[1][0] = 0
        B[2][0] = (c * 0.25 * (s - 1) - d * 0.25 * (t - 1)) / DetJacobian(s, t)

        B[0][1] = 0
        B[1][1] = (c * 0.25 * (s - 1) - d * 0.25 * (t - 1)) / DetJacobian(s, t)
        B[2][1] = (a * 0.25 * (t - 1) - b * 0.25 * (s - 1)) / DetJacobian(s, t)

        B[0][2] = (a * 0.25 * (1 - t) - b * 0.25 * (-1 - s)) / DetJacobian(s, t)
        B[1][2] = 0
        B[2][2] = (c * 0.25 * (-1 - s) - d * 0.25 * (1 - t)) / DetJacobian(s, t)

        B[0][3] = 0
        B[1][3] = (c * 0.25 * (-1 - s) - d * 0.25 * (1 - t)) / DetJacobian(s, t)
        B[2][3] = (a * 0.25 * (1 - t) - b * 0.25 * (-1 - s)) / DetJacobian(s, t)

        B[0][4] = (a * 0.25 * (1 + t) - b * 0.25 * (1 + s)) / DetJacobian(s, t)
        B[1][4] = 0
        B[2][4] = (c * 0.25 * (1 + s) - d * 0.25 * (1 + t)) / DetJacobian(s, t)

        B[0][5] = 0
        B[1][5] = (c * 0.25 * (1 + s) - d * 0.25 * (1 + t)) / DetJacobian(s, t)
        B[2][5] = (a * 0.25 * (1 + t) - b * 0.25 * (1 + s)) / DetJacobian(s, t)

        B[0][6] = (a * 0.25 * (-1 - t) - b * 0.25 * (1 - s)) / DetJacobian(s, t)
        B[1][6] = 0
        B[2][6] = (c * 0.25 * (1 - s) - d * 0.25 * (-1 - t)) / DetJacobian(s, t)

        B[0][7] = 0
        B[1][7] = (c * 0.25 * (1 - s) - d * 0.25 * (-1 - t)) / DetJacobian(s, t)
        B[2][7] = (a * 0.25 * (-1 - t) - b * 0.25 * (1 - s)) / DetJacobian(s, t)

        return B


    Ke1 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(s1, t1)), matrixmult(D, BMatrix(s1, t1))), DetJacobian(s1, t1)),
        Thickness)
    Ke2 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(s2, t2)), matrixmult(D, BMatrix(s2, t2))), DetJacobian(s2, t2)),
        Thickness)
    Ke3 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(s3, t3)), matrixmult(D, BMatrix(s3, t3))), DetJacobian(s3, t3)),
        Thickness)
    Ke4 = matrixscalarmult(
        matrixscalarmult(matrixmult(transpose(BMatrix(s4, t4)), matrixmult(D, BMatrix(s4, t4))), DetJacobian(s4, t4)),
        Thickness)

    Ke = matrixsum(matrixsum(matrixsum(Ke1, Ke2), Ke3), Ke4)


    # ------------------------------------- MASS Matrix -------------------------------------------------
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


    MassMatrix1 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(s1, t1)), NMatrix(s1, t1)), DetJacobian(s1, t1)),
                         Thickness), density)
    MassMatrix2 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(s2, t2)), NMatrix(s2, t2)), DetJacobian(s2, t2)),
                         Thickness), density)
    MassMatrix3 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(s3, t3)), NMatrix(s3, t3)), DetJacobian(s3, t3)),
                         Thickness), density)
    MassMatrix4 = matrixscalarmult(
        matrixscalarmult(matrixscalarmult(matrixmult(transpose(NMatrix(s4, t4)), NMatrix(s4, t4)), DetJacobian(s4, t4)),
                         Thickness), density)

    MassMatrixelement = matrixsum(matrixsum(matrixsum(MassMatrix1, MassMatrix2), MassMatrix3), MassMatrix4)

    # *************************************************************************************
    # GLOBAL STIFFNESS ANS MASS MATRIX ASSEMBLY
    # *************************************************************************************
    for r in range(0, 4):
        for t in range(0, 2):
            for z in range(0, 4):
                for m in range(0, 2):
                    rw = NDOFN * (EL_Data[i - 1][r + 1] - 1) + t
                    cl = NDOFN * (EL_Data[i - 1][z + 1] - 1) + m
                    KG[rw][cl] = KG[rw][cl] + Ke[NDOFN * r + t][NDOFN * z + m]
                    Mmatrix[rw][cl] = Mmatrix[rw][cl] + MassMatrixelement[NDOFN * r + t][NDOFN * z + m]

print('GLOBAL STIFFNESS MATRIX READY')

# *************************************************************************************
# BODY FORCE + LOADS
# *************************************************************************************
# LOADS_total=np.add(matrixmult(Mmatrix,Acceleration),LOADS)
LOADS_total = LOADS
# --------------removing the body forces for the nodes with BCs
index = 0
while index < (NDOF):
    if BC[index][0] == 1:
        LOADS_total[index][0] = 0
    else:
        LOADS_total[index][0] = LOADS_total[index][0]
    index = index + 1
# print(LOADS_total)

# *************************************************************************************
# APPLYING BOUNDARY CONDITIONS
# *************************************************************************************
# --------------Filling Global AUXILIAR Stiffness Matrix with zeros
KG2 = [[0 for row in range(NDOF)] for col in range(NDOF)]
# ----------------------------------------------------------------------------
index = 0
while index < (NDOF):
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
print('BOUNDARY CONDITIONS APPLIED')
# *************************************************************************************
# Solving {U} = [K]^(-1) x {F}
# *************************************************************************************
# ------------------------- Inverting the KG2 Matrix ---------------------
print('INVERTING GLOBAL STIFFNESS MATRIX')
KGINV = np.linalg.pinv(KG2)  # --------------- used Moore-Penrose pseudo-inverse scheme
# KGINV=invert(KG2)
# print (KGINV)
print('GLOBAL STIFFNESS MATRIX INVERTED')

# ------------------------- Calculating the nodal displacements (U = [K]^(-1) x {F}) ---------------------
Displacement = matrixmult(KGINV, LOADS_total)



#maxm nodal displacement
d_fem = max(Displacement[:])[0]
print(d_fem)

def savePlots(node_list, U):
    x = []
    y = []
    u = []
    v = []
    for i in range(0, len(node_list)):
        x.append(node_list[i][1])
        y.append(node_list[i][2])
        u.append(U[2 * i][0])
        v.append(U[2 * i + 1][0])

    plt.figure(figsize=(12, 8))
    plt.quiver(x, y, u, v)
    plt.show()
    plt.clf()


savePlots(Node_Data, Displacement)


