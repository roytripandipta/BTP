
import matplotlib.pyplot as plt
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
def create_nodes(length, breadth, nodes_in_x, nodes_in_y):
    x = np.linspace(0, length, nodes_in_x)
    y = np.linspace(0, breadth, nodes_in_y)
    xx, yy = np.meshgrid(x, y)
    # dist = x[1]-x[0]
    mid = len(xx[0]) // 2
    crack_l = nodes_in_y // 4
    x1 = xx[:crack_l]
    x2 = xx[crack_l:].tolist()
    mid = len(xx[0]) // 2
    x1 = np.insert(x1, mid + 1, x1[:, mid], axis=1).tolist()
    lst = x2 + x1
    y1 = yy[:-crack_l].tolist()
    y2 = yy[-crack_l:]
    mid = len(yy[0]) // 2
    y2 = np.insert(y2, mid + 1, y2[:, mid], axis=1).tolist()
    lst1 = y1 + y2
    no_of_nodes = 0
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            no_of_nodes += 1

    arr = np.zeros((no_of_nodes + 1, 3))
    count = 1
    for j in range(len(lst)):
        for k in range(len(lst[j])):
            # arr[count][2] = yy[j][k]
            arr[count][0] = int(count)
            arr[count][1] = lst[j][k]
            arr[count][2] = lst1[j][k]
            count += 1
    arr = arr[1:].tolist()
    arr = [[int(x[0]), x[1], x[2]] for x in arr]
    return arr, no_of_nodes, lst, lst1


def create_connectivity(length, breadth, nodes_in_x, nodes_in_y):
    _, _, l1, l2 = create_nodes(length, breadth, nodes_in_x, nodes_in_y)
    g = []
    count = 1
    for i in range(len(l1)):
        lt = []
        for j in range(len(l1[i])):
            lt.append(count)
            count += 1
        g.append(lt)
    eles = []
    cnt = 1
    for i in range(len(g) - 1):
        flag = False
        for j in range(len(g[i]) - 1):
            l = []
            if l1[i][j] != l1[i][j + 1]:
                if (l1[i + 1][j] == l1[i + 1][j + 1]) or flag:
                    l.append(cnt)
                    l.append(g[i][j])
                    l.append(g[i][j + 1])
                    l.append(g[i + 1][j + 2])
                    l.append(g[i + 1][j + 1])
                    eles.append(l)
                    flag = True
                else:
                    l.append(cnt)
                    l.append(g[i][j])
                    l.append(g[i][j + 1])
                    l.append(g[i + 1][j + 1])
                    l.append(g[i + 1][j])
                    eles.append(l)
                cnt += 1
    return eles, len(eles), g


def extract_elem(g, eles, height, arr, disp):
    num = 0
    for i in range(len(g)):
        if len(g[i]) < len(g[i + 1]):
            pos = len(g[i]) // 2
            num = g[i][pos]
            break
    element = []
    for i in range(len(eles)):
        if num in eles[i][1:]:
            element.append(eles[i])
    crack_points = []
    is_present = []
    arr1 = deepcopy(arr)
    for i in range(len(element)):
        for j in element[i][1:]:
            if j not in is_present:
                arr1[j - 1].extend(
                    [0.0, arr1[j - 1][1] + disp[2 * (j - 1)][0], arr1[j - 1][2] + disp[2 * (j - 1) + 1][0], 0.0])
                crack_points.append(arr1[j - 1])
                is_present.append(j)
    crack = []
    crack.extend([arr[num - 1][1], float(height)])
    crack.extend(arr[num - 1][1:])
    return crack, element, crack_points


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
def apply_loads(Node_Data, length, F):
    NNodes = len(Node_Data)
    NDOF = 2 * NNodes
    l = 0
    for x_ in Node_Data:
        if x_[1] == 0:
            l += 1
    F = F / l
    LOADS = [[0 for row in range(1)] for col in range(NDOF)]  # required matrix for loading conditions
    for x in Node_Data:
        if x[1] == length:
            LOADS[2 * (x[0] - 1)][0] = F
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


def theoretical_disp(length, breadth, t, E, P):
    A = breadth * t
    d_th = (P * length) / (A * E)
    return d_th


def showPlots(node_list, U):
    x = []
    y = []
    u = []
    v = []
    for j in range(0, len(node_list)):
        x.append(node_list[j][1])
        y.append(node_list[j][2])
        u.append(U[2 * j][0])
        v.append(U[2 * j + 1][0])
    plt.figure()
    plt.quiver(x, y, u, v)
    plt.show()
    plt.close()


# post processing
def compute_db(D, ls, lt, X, Y):
    DB1 = matrixmult(D, BMatrix(ls[0], lt[0], X, Y))
    DB2 = matrixmult(D, BMatrix(ls[1], lt[1], X, Y))
    DB3 = matrixmult(D, BMatrix(ls[2], lt[2], X, Y))
    DB4 = matrixmult(D, BMatrix(ls[3], lt[3], X, Y))

    return DB1, DB2, DB3, DB4


def compute_dispElement(Displacement, EL_Data, i):
    Disp_elem = [[0 for row in range(1)] for col in range(8)]

    Disp_elem[0][0] = Displacement[2 * (EL_Data[i - 1][1]) - 2][0]
    Disp_elem[1][0] = Displacement[2 * (EL_Data[i - 1][1]) - 1][0]
    Disp_elem[2][0] = Displacement[2 * (EL_Data[i - 1][2]) - 2][0]
    Disp_elem[3][0] = Displacement[2 * (EL_Data[i - 1][2]) - 1][0]
    Disp_elem[4][0] = Displacement[2 * (EL_Data[i - 1][3]) - 2][0]
    Disp_elem[5][0] = Displacement[2 * (EL_Data[i - 1][3]) - 1][0]
    Disp_elem[6][0] = Displacement[2 * (EL_Data[i - 1][4]) - 2][0]
    Disp_elem[7][0] = Displacement[2 * (EL_Data[i - 1][4]) - 1][0]
    return Disp_elem


def compute_strain_integration_point(ls, lt, X, Y, Disp_elem):
    Strains_integrationpoint1 = matrixmult(BMatrix(ls[0], lt[0], X, Y), Disp_elem)
    Strains_integrationpoint2 = matrixmult(BMatrix(ls[1], lt[1], X, Y), Disp_elem)
    Strains_integrationpoint3 = matrixmult(BMatrix(ls[2], lt[2], X, Y), Disp_elem)
    Strains_integrationpoint4 = matrixmult(BMatrix(ls[3], lt[3], X, Y), Disp_elem)
    return Strains_integrationpoint1, Strains_integrationpoint2, Strains_integrationpoint3, Strains_integrationpoint4


def compute_stress_integration_point(db, Disp_elem):
    Stress_integrationpoint1 = matrixmult(db[0], Disp_elem)
    Stress_integrationpoint2 = matrixmult(db[1], Disp_elem)
    Stress_integrationpoint3 = matrixmult(db[2], Disp_elem)
    Stress_integrationpoint4 = matrixmult(db[3], Disp_elem)
    return Stress_integrationpoint1, Stress_integrationpoint2, Stress_integrationpoint3, Stress_integrationpoint4


def compute_Sx(Stress_integrationPoint):
    StressX = [[0 for row in range(1)] for col in range(4)]

    StressX[0][0] = Stress_integrationPoint[0][0][0]
    StressX[1][0] = Stress_integrationPoint[1][0][0]
    StressX[2][0] = Stress_integrationPoint[2][0][0]
    StressX[3][0] = Stress_integrationPoint[3][0][0]

    return StressX


def compute_Sy(Stress_integrationPoint):
    StressY = [[0 for row in range(1)] for col in range(4)]

    StressY[0][0] = Stress_integrationPoint[0][1][0]
    StressY[1][0] = Stress_integrationPoint[1][1][0]
    StressY[2][0] = Stress_integrationPoint[2][1][0]
    StressY[3][0] = Stress_integrationPoint[3][1][0]

    return StressY


def compute_Sxy(Stress_integrationPoint):
    StressXY = [[0 for row in range(1)] for col in range(4)]

    StressXY[0][0] = Stress_integrationPoint[0][2][0]
    StressXY[1][0] = Stress_integrationPoint[1][2][0]
    StressXY[2][0] = Stress_integrationPoint[2][2][0]
    StressXY[3][0] = Stress_integrationPoint[3][2][0]

    return StressXY


def compute_principal1(StressX, StressY, StressXY):
    StressPrinc1 = [[0 for row in range(1)] for col in range(4)]

    StressPrinc1[0][0] = (StressX[0][0] + StressY[0][0]) / 2 + (
            ((StressX[0][0] - StressY[0][0]) / 2) ** 2 + (StressXY[0][0]) ** 2) ** 0.5
    StressPrinc1[1][0] = (StressX[1][0] + StressY[1][0]) / 2 + (
            ((StressX[1][0] - StressY[1][0]) / 2) ** 2 + (StressXY[1][0]) ** 2) ** 0.5
    StressPrinc1[2][0] = (StressX[2][0] + StressY[2][0]) / 2 + (
            ((StressX[2][0] - StressY[2][0]) / 2) ** 2 + (StressXY[2][0]) ** 2) ** 0.5
    StressPrinc1[3][0] = (StressX[3][0] + StressY[3][0]) / 2 + (
            ((StressX[3][0] - StressY[3][0]) / 2) ** 2 + (StressXY[3][0]) ** 2) ** 0.5

    return StressPrinc1


def compute_principal2(StressX, StressY, StressXY):
    StressPrinc2 = [[0 for row in range(1)] for col in range(4)]
    StressPrinc2[0][0] = (StressX[0][0] + StressY[0][0]) / 2 - (
            ((StressX[0][0] - StressY[0][0]) / 2) ** 2 + (StressXY[0][0]) ** 2) ** 0.5
    StressPrinc2[1][0] = (StressX[1][0] + StressY[1][0]) / 2 - (
            ((StressX[1][0] - StressY[1][0]) / 2) ** 2 + (StressXY[1][0]) ** 2) ** 0.5
    StressPrinc2[2][0] = (StressX[2][0] + StressY[2][0]) / 2 - (
            ((StressX[2][0] - StressY[2][0]) / 2) ** 2 + (StressXY[2][0]) ** 2) ** 0.5
    StressPrinc2[3][0] = (StressX[3][0] + StressY[3][0]) / 2 - (
            ((StressX[3][0] - StressY[3][0]) / 2) ** 2 + (StressXY[3][0]) ** 2) ** 0.5

    return StressPrinc2


def compute_VM(StressPrinc1, StressPrinc2):
    StressVM = [[0 for row in range(1)] for col in range(4)]

    StressVM[0][0] = ((StressPrinc1[0][0]) ** 2 - StressPrinc1[0][0] * StressPrinc2[0][0] + (
        StressPrinc2[0][0]) ** 2) ** 0.5
    StressVM[1][0] = ((StressPrinc1[1][0]) ** 2 - StressPrinc1[1][0] * StressPrinc2[1][0] + (
        StressPrinc2[1][0]) ** 2) ** 0.5
    StressVM[2][0] = ((StressPrinc1[2][0]) ** 2 - StressPrinc1[2][0] * StressPrinc2[2][0] + (
        StressPrinc2[2][0]) ** 2) ** 0.5
    StressVM[3][0] = ((StressPrinc1[3][0]) ** 2 - StressPrinc1[3][0] * StressPrinc2[3][0] + (
        StressPrinc2[3][0]) ** 2) ** 0.5

    return StressVM


# def extrapolate_result(EXTRAP, StressX, StressY, StressXY, StressVM, StressPrinc1, StressPrinc2):
#     StressX_EXTRAP = matrixmult(EXTRAP, StressX)
#     StressY_EXTRAP = matrixmult(EXTRAP, StressY)
#     StressXY_EXTRAP = matrixmult(EXTRAP, StressXY)
#     StressPrinc1_EXTRAP = matrixmult(EXTRAP, StressPrinc1)
#     StressPrinc2_EXTRAP = matrixmult(EXTRAP, StressPrinc2)
#     StressVM_EXTRAP = matrixmult(EXTRAP, StressVM)
#
#     return StressX_EXTRAP, StressY_EXTRAP, StressXY_EXTRAP, StressPrinc1_EXTRAP, StressPrinc2_EXTRAP, StressVM_EXTRAP

"""
def paraview_result(i, Node_Data, EL_Data, StressX, StressY, StressXY, StressPrinc1, StressPrinc2, StressVM,
                    StressX_EXTRAP, StressY_EXTRAP, StressXY_EXTRAP, StressPrinc1_EXTRAP, StressPrinc2_EXTRAP,
                    StressVM_EXTRAP, Strains_integrationpoint):
    NNodes = 2 * Node_Data
    PARAVIEW_RESULTS = [[None for row in range(72)] for col in range(NNodes)]
    for r in range(0, NNodes):
        PARAVIEW_RESULTS[r][0] = Node_Data[r][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][1] = StressX[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][2] = StressX[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][3] = StressX[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][4] = StressX[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][6] = StressY[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][7] = StressY[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][8] = StressY[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][9] = StressY[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][11] = StressXY[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][12] = StressXY[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][13] = StressXY[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][14] = StressXY[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][16] = StressPrinc1[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][17] = StressPrinc1[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][18] = StressPrinc1[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][19] = StressPrinc1[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][21] = StressPrinc2[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][22] = StressPrinc2[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][23] = StressPrinc2[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][24] = StressPrinc2[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][26] = StressVM[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][27] = StressVM[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][28] = StressVM[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][29] = StressVM[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][36] = StressX_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][37] = StressX_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][38] = StressX_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][39] = StressX_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][41] = StressY_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][42] = StressY_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][43] = StressY_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][44] = StressY_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][46] = StressXY_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][47] = StressXY_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][48] = StressXY_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][49] = StressXY_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][51] = StressPrinc1_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][52] = StressPrinc1_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][53] = StressPrinc1_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][54] = StressPrinc1_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][56] = StressPrinc2_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][57] = StressPrinc2_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][58] = StressPrinc2_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][59] = StressPrinc2_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][61] = StressVM_EXTRAP[0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][62] = StressVM_EXTRAP[1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][63] = StressVM_EXTRAP[2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][64] = StressVM_EXTRAP[3][0]

    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][66] = Strains_integrationpoint[0][0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][66] = Strains_integrationpoint[1][0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][66] = Strains_integrationpoint[2][0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][66] = Strains_integrationpoint[3][0][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][70] = Strains_integrationpoint[0][1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][70] = Strains_integrationpoint[1][1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][70] = Strains_integrationpoint[2][1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][70] = Strains_integrationpoint[3][1][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][1] - 1][71] = Strains_integrationpoint[0][2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][2] - 1][71] = Strains_integrationpoint[1][2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][3] - 1][71] = Strains_integrationpoint[2][2][0]
    PARAVIEW_RESULTS[EL_Data[i - 1][4] - 1][71] = Strains_integrationpoint[3][2][0]

    # ================================= averaging the results ====================

    return PARAVIEW_RESULTS

"""


def make_average(PARAVIEW_RESULTS, NNodes):
    print('PARAVIEW FILE OUTPUT')
    # ------------------------------------ SX
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(1, 5):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][5] = Sum / k
    # ------------------------------------ SY
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(6, 10):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][10] = Sum / k
    # ------------------------------------ SXY
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(11, 15):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][15] = Sum / k
    # ------------------------------------ S1
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(16, 20):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][20] = Sum / k

    # ------------------------------------ S2
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(21, 25):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][25] = Sum / k

    # ------------------------------------ VM
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(26, 30):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][30] = Sum / k
    # ------------------------------------ SX EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(36, 40):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][40] = Sum / k

    # ------------------------------------ SY EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(41, 45):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][45] = Sum / k

    # ------------------------------------ SXY EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(46, 50):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][50] = Sum / k
    # ------------------------------------ S1 EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(51, 55):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][55] = Sum / k
    # ------------------------------------ S2 EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(56, 60):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][60] = Sum / k
    # ------------------------------------ VM EXTRAPOLATED
    for i in range(0, NNodes):
        k = 0
        Sum = 0
        for j in range(61, 65):
            if PARAVIEW_RESULTS[i][j] is not None:
                k = k + 1
                Sum = Sum + PARAVIEW_RESULTS[i][j]
        if k == 0:
            print(i)
        PARAVIEW_RESULTS[i][65] = Sum / k
    # ------------------------------------ VECTORS RESULTS
    for i in range(0, NNodes):
        PARAVIEW_RESULTS[i][31] = np.arctan(
            (2 * PARAVIEW_RESULTS[i][15]) / (PARAVIEW_RESULTS[i][5] - PARAVIEW_RESULTS[i][10])) / 2 * 180 / np.pi
    for i in range(0, NNodes):
        PARAVIEW_RESULTS[i][32] = np.cos(PARAVIEW_RESULTS[i][31] / (180 / np.pi)) * PARAVIEW_RESULTS[i][20]
    for i in range(0, NNodes):
        PARAVIEW_RESULTS[i][33] = np.sin(PARAVIEW_RESULTS[i][31] / (180 / np.pi)) * PARAVIEW_RESULTS[i][20]
    for i in range(0, NNodes):
        PARAVIEW_RESULTS[i][34] = np.sin(PARAVIEW_RESULTS[i][31] / (180 / np.pi)) * PARAVIEW_RESULTS[i][25]
    for i in range(0, NNodes):
        PARAVIEW_RESULTS[i][35] = np.cos(PARAVIEW_RESULTS[i][31] / (180 / np.pi)) * PARAVIEW_RESULTS[i][25]

    return PARAVIEW_RESULTS
