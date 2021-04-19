import fem_utilities as fem
import numpy as np

Analysis_Type = 1
E = 210e9
Poisson = 0
Thickness = 1
Material_data = [[Analysis_Type, E, Poisson, Thickness]]
pi = np.arccos(-1)
if Analysis_Type == 1:
    print('Type of Analysis: Plane Stress')
else:
    print('Type of Analysis: Plane Strain')

Node_Data = fem.make_nodes(10, 20, 3, 5)
print(Node_Data)
EL_Data, NElements = fem.connectivity(10, 20, 3, 5)
print(EL_Data)
print(NElements)
BC = fem.make_BC(Node_Data)
LOADS = fem.apply_loads(Node_Data, 10, 10000)
NDOFN = 2
NNodes = len(Node_Data)
NDOF = 2 * NNodes
density = 7800

# ----------------- Creating an array with node data

# print(Acceleration)
# ======================================================================
#
#             READING ALL DATA FROM EXTERNAL FILE - END
#
# ======================================================================



# ===================================================================================
#
#                            FUNCTIONS DEFINITION
#
# ===================================================================================

# --------------------------------------- MATRIX MULTIPLICATION FUNCTION ---------------
def matrixmult(A, B):
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    return C


# --------------------------------------- MATRIX SUM FUNCTION ---------------
def matrixsum(A, B):
    C = [[0 for row in range(len(A))] for col in range(len(A[0]))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            C[i][j] = A[i][j] + B[i][j]
    return C


# --------------------------------------- MATRIX MULTIPLICATION FUNCTION by an SCALAR ---------------
def matrixscalarmult(A, b):
    for row in A:
        Multiscalar = [[A[i][j] * b for i in range(len(A))] for j in range(len(A[0]))]
    return Multiscalar


# --------------------------------------- MATRIX TRANSPOSE FUNCTION ---------------
def transpose(A):
    for row in A:
        Trans = [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]
    return Trans


# ===================================================================================
#
#
#
#                                      FEA CODE
#
#
#
# ===================================================================================
# ----------------------------Integration Points Coordinates (4 integration points)
s1 = -0.577350269
t1 = -0.577350269
s2 = 0.577350269
t2 = -0.577350269
s3 = 0.577350269
t3 = 0.577350269
s4 = -0.577350269
t4 = 0.577350269
# -----------------------------------Filling Global Stiffness Matrix with zeros
KG = [[0 for row in range(NDOF)] for col in range(NDOF)]
LOADS_total = [[0 for row in range(1)] for col in range(NDOF)]
# -----------------------------------Filling Global Mass Matrix with zeros
Mmatrix = [[0 for row in range(NDOF)] for col in range(NDOF)]
# -----------------------------------Filling Results Matrix for each element with zeros
Disp_elem = [[0 for row in range(1)] for col in range(8)]
StressX = [[0 for row in range(1)] for col in range(4)]
StressY = [[0 for row in range(1)] for col in range(4)]
StressXY = [[0 for row in range(1)] for col in range(4)]
StressPrinc1 = [[0 for row in range(1)] for col in range(4)]
StressPrinc2 = [[0 for row in range(1)] for col in range(4)]
StressVM = [[0 for row in range(1)] for col in range(4)]
StressX_EXTRAP = [[0 for row in range(1)] for col in range(4)]
StressY_EXTRAP = [[0 for row in range(1)] for col in range(4)]
StressXY_EXTRAP = [[0 for row in range(1)] for col in range(4)]
StressPrinc1_EXTRAP = [[0 for row in range(1)] for col in range(4)]
StressPrinc2_EXTRAP = [[0 for row in range(1)] for col in range(4)]
StressVM_EXTRAP = [[0 for row in range(1)] for col in range(4)]
Strains_integrationpoint1 = [[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint2 = [[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint3 = [[0 for row in range(1)] for col in range(3)]
Strains_integrationpoint4 = [[0 for row in range(1)] for col in range(3)]
PARAVIEW_RESULTS = [[None for row in range(72)] for col in range(NNodes)]
for r in range(0, NNodes):
    PARAVIEW_RESULTS[r][0] = Node_Data[r][0]
# ------------------- nodal extrapolation matrix -------------------------
EXTRAP = [[1 + 0.5 * 3 ** 0.5, -0.5, 1 - 0.5 * 3 ** 0.5, -0.5], [-0.5, 1 + 0.5 * 3 ** 0.5, -0.5, 1 - 0.5 * 3 ** 0.5],
          [1 - 0.5 * 3 ** 0.5, -0.5, 1 + 0.5 * 3 ** 0.5, -0.5], [-0.5, 1 - 0.5 * 3 ** 0.5, -0.5, 1 + 0.5 * 3 ** 0.5]]
# ==================================Defining material behavior - Plane stress or Plane strain
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
                    KG[rw][cl] = KG[rw][cl] + Ke[NDOFN * (r) + t][NDOFN * (z) + m]
                    Mmatrix[rw][cl] = Mmatrix[rw][cl] + MassMatrixelement[NDOFN * (r) + t][NDOFN * (z) + m]

# print (KG)
print('GLOBAL STIFFNESS MATRIX READY')
# print(Mmatrix)
# *************************************************************************************
# BODY FORCE + LOADS
# *************************************************************************************
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

print(Displacement)
