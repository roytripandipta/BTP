import fem_utilities as fem
import numpy as np
from tkinter import *


# Note all units are in SI
def compute_stiffness(Node_Data, EL_Data, NElements, ls, lt, D, Thickness, density, NDOFN):
    NDOF = 2 * len(Node_Data)
    Mmatrix = [[0 for row in range(NDOF)] for col in range(NDOF)]
    KG = [[0 for row in range(NDOF)] for col in range(NDOF)]
    for i in range(1, NElements + 1):
        X = [[Node_Data[EL_Data[i - 1][1] - 1][1], Node_Data[EL_Data[i - 1][2] - 1][1],
              Node_Data[EL_Data[i - 1][3] - 1][1], Node_Data[EL_Data[i - 1][4] - 1][1]]]
        Y = [[Node_Data[EL_Data[i - 1][1] - 1][2], Node_Data[EL_Data[i - 1][2] - 1][2],
              Node_Data[EL_Data[i - 1][3] - 1][2], Node_Data[EL_Data[i - 1][4] - 1][2]]]
        Ke = fem.make_Ke(ls, lt, X, Y, D, Thickness)

        MassMatrixelement = fem.mass_matrix(ls, lt, X, Y, Thickness, density)
        for r in range(0, 4):
            for t in range(0, 2):
                for z in range(0, 4):
                    for m in range(0, 2):
                        rw = NDOFN * (EL_Data[i - 1][r + 1] - 1) + t
                        cl = NDOFN * (EL_Data[i - 1][z + 1] - 1) + m
                        KG[rw][cl] = KG[rw][cl] + Ke[NDOFN * r + t][NDOFN * z + m]
                        Mmatrix[rw][cl] = Mmatrix[rw][cl] + MassMatrixelement[NDOFN * r + t][NDOFN * z + m]

    return KG, Mmatrix


# print(disp)
#
# print(len(KG2[0]))


def solve():
    Analysis_Type = 1  # 1 - plain stress, 2 - plain strain
    E = 210e9  # Young modulus
    Poisson = 0.3
    density = 7800
    NDOFN = 2
    Thickness = 1
    AccelX = 0
    AccelY = 0
    length = length.get()
    height = height.get()
    no_of_nodes_in_x = 3
    no_of_nodes_in_y = 5
    # make a 2d matrix for material data
    Material_data = [[Analysis_Type, E, Poisson, Thickness]]

    pi = np.arccos(-1)
    if Analysis_Type == 1:
        print('Type of Analysis: Plane Stress')
    else:
        print('Type of Analysis: Plane Strain')
    ls = [-0.577350269, 0.577350269, 0.577350269, -0.577350269]
    lt = [-0.577350269, -0.577350269, 0.577350269, 0.577350269]

    try:
        Node_Data = fem.make_nodes(length, breadth, no_of_nodes_in_x, no_of_nodes_in_y)
        # print(Node_Data)
        EL_Data, NElements = fem.connectivity(length, breadth, no_of_nodes_in_x, no_of_nodes_in_y)
        # print(EL_Data)
        BC = fem.make_BC(Node_Data)
        # print(BC)
        LOADS = fem.apply_loads(Node_Data, length)

        # print(LOADS)
        D = fem.calculate_D(Analysis_Type, E, Poisson)
        # print(D)
        NDOF = 2 * len(Node_Data)
        KG, _ = compute_stiffness(Node_Data, EL_Data, NElements, ls, lt, D, Thickness, density, 2)
        LOADS_total = fem.remove_force_from_BC(LOADS, BC, NDOF)
        KG2 = fem.apply_BC(BC, KG, NDOF)
        disp = fem.compute_displacement(KG2, LOADS_total)
        print(fem.max_disp(disp))
        print(fem.calc_average(disp, Node_Data, length))
        return disp
    except Exception as e:
        print(f"error in solution : {e}")

        return None


if __name__ == "__main__":
    root = Tk()
    root.geometry("644x344")
    Label(root, text="Welcome to FEM World", font="comicsansms 13 bold", pady=15).grid(row=0, column=4)
    # Text for our form
    disp = None
    Label(root, text="Length").grid(row=1, column=2)
    Label(root, text="Height").grid(row=2, column=2)
    Label(root, text="Thickness").grid(row=3, column=2)
    Label(root, text="Poisson Ratio").grid(row=4, column=2)
    Label(root, text="Young Modulus").grid(row=5, column=2)
    Label(root, text="No of nodes in x direction").grid(row=6, column=2)
    Label(root, text="No of nodes in y direction").grid(row=7, column=2)

    length = DoubleVar()
    height = DoubleVar()
    thickness = DoubleVar()
    poisson_ratio = DoubleVar()
    E = DoubleVar()
    nodes_x = IntVar()
    nodes_y = IntVar()

    Entry(root, textvariable=length).grid(row=1, column=4)
    Entry(root, textvariable=height).grid(row=2, column=4)
    Entry(root, textvariable=thickness).grid(row=3, column=4)
    Entry(root, textvariable=poisson_ratio).grid(row=4, column=4)
    Entry(root, textvariable=E).grid(row=5, column=4)
    Entry(root, textvariable=nodes_x).grid(row=6, column=4)
    Entry(root, textvariable=nodes_y).grid(row=7, column=4)

    Submit = Button(text="Submit", command = lambda:disp == solve()).grid(row=8, column=4, pady=15)
    root.mainloop()
