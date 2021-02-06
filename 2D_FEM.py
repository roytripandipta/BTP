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
# Text for our form
def run():
    root = Tk()
    root.geometry("500x500")
    root.minsize(500, 500)
    Label(root, text="Welcome to FEM World", font="comicsansms 13 bold", pady=15).grid(row=0, column=4)
    root.title("FEM Calculator")
    root.wm_iconbitmap("icon.ico")

    def solve():
        Analysis_Type = 1  # 1 - plain stress, 2 - plain strainq
        E = Young_mod.get()  # Young modulus
        Poisson = Poisson_ratio.get()
        density = 7800
        NDOFN = 2
        Thickness = Thickness_.get()
        AccelX = 0
        AccelY = 0
        length = Length.get()
        breadth = Height.get()
        F = Force.get()
        no_of_nodes_in_x = nodes_x.get()
        no_of_nodes_in_y = nodes_y.get()
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
            LOADS = fem.apply_loads(Node_Data, length, F)

            # print(LOADS)
            D = fem.calculate_D(Analysis_Type, E, Poisson)
            # print(D)
            NDOF = 2 * len(Node_Data)
            KG, _ = compute_stiffness(Node_Data, EL_Data, NElements, ls, lt, D, Thickness, density, 2)
            LOADS_total = fem.remove_force_from_BC(LOADS, BC, NDOF)
            KG2 = fem.apply_BC(BC, KG, NDOF)
            displacement = fem.compute_displacement(KG2, LOADS_total)
            avg_d = fem.calc_average(displacement, Node_Data, length)
            print(avg_d)
            scvalue1.set(avg_d)
            # print(displacement)
            max_d = fem.max_disp(displacement)
            print(max_d)
            scvalue2.set(max_d)
            # screen2.update()
            d_th = fem.theoretical_disp(length, breadth, Thickness, E, F)
            scvalue3.set(d_th)
            fem.showPlots(Node_Data, displacement)
            np.savetxt("result/displacement.csv", displacement)
            NNodes = len(Node_Data)
            PARAVIEW_RESULTS = [[None for row in range(72)] for col in range(NNodes)]
            for i in range(1, NElements + 1):
                X = [[Node_Data[EL_Data[i - 1][1] - 1][1], Node_Data[EL_Data[i - 1][2] - 1][1],
                      Node_Data[EL_Data[i - 1][3] - 1][1], Node_Data[EL_Data[i - 1][4] - 1][1]]]
                Y = [[Node_Data[EL_Data[i - 1][1] - 1][2], Node_Data[EL_Data[i - 1][2] - 1][2],
                      Node_Data[EL_Data[i - 1][3] - 1][2], Node_Data[EL_Data[i - 1][4] - 1][2]]]

                D = fem.calculate_D(Analysis_Type, E, Poisson)
                DB = fem.compute_db(D, ls, lt, X, Y)
                Disp_elem = fem.compute_dispElement(displacement, EL_Data, i)
                Strains_integrationpoint = fem.compute_strain_integration_point(ls, lt, X, Y, Disp_elem)
                Stress_integrationpoint = fem.compute_stress_integration_point(DB, Disp_elem)
                StressX = fem.compute_Sx(Stress_integrationpoint)
                # print(StressX)
                StressY = fem.compute_Sy(Stress_integrationpoint)
                StressXY = fem.compute_Sxy(Stress_integrationpoint)
                StressPrinc1 = fem.compute_principal1(StressX, StressY, StressXY)
                StressPrinc2 = fem.compute_principal2(StressX, StressY, StressXY)
                StressVM = fem.compute_VM(StressPrinc1, StressPrinc2)
                # print(Strains_integrationpoint)
                EXTRAP = [[1 + 0.5 * 3 ** 0.5, -0.5, 1 - 0.5 * 3 ** 0.5, -0.5],
                          [-0.5, 1 + 0.5 * 3 ** 0.5, -0.5, 1 - 0.5 * 3 ** 0.5],
                          [1 - 0.5 * 3 ** 0.5, -0.5, 1 + 0.5 * 3 ** 0.5, -0.5],
                          [-0.5, 1 - 0.5 * 3 ** 0.5, -0.5, 1 + 0.5 * 3 ** 0.5]]
                StressX_EXTRAP = fem.matrixmult(EXTRAP, StressX)
                StressY_EXTRAP = fem.matrixmult(EXTRAP, StressY)
                StressXY_EXTRAP = fem.matrixmult(EXTRAP, StressXY)
                StressPrinc1_EXTRAP = fem.matrixmult(EXTRAP, StressPrinc1)
                StressPrinc2_EXTRAP = fem.matrixmult(EXTRAP, StressPrinc2)
                StressVM_EXTRAP = fem.matrixmult(EXTRAP, StressVM)
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
            P = fem.make_average(PARAVIEW_RESULTS, NNodes)
            s = 0
            for i in range(0, NNodes):
                s += P[i][55]
                print(P[i][55])
            print(s / NNodes)

        except Exception as e:
            print(f"error in solution : {e}")

    Label(root, text="Length").grid(row=1, column=2)
    Label(root, text="Height").grid(row=2, column=2)
    Label(root, text="Thickness").grid(row=3, column=2)
    Label(root, text="Poisson Ratio").grid(row=4, column=2)
    Label(root, text="Young Modulus").grid(row=5, column=2)
    Label(root, text="Force").grid(row=6, column=2)
    Label(root, text="No of nodes in x direction").grid(row=7, column=2)
    Label(root, text="No of nodes in y direction").grid(row=8, column=2)
    Label(root, text="Average Displacement", font="lucida 10 bold").grid(row=10, column=2)
    Label(root, text="Maximum Displacement", font="lucida 10 bold").grid(row=11, column=2)
    Label(root, text="Theoretical Displacement", font="lucida 10 bold").grid(row=12, column=2)
    Length = DoubleVar()
    Height = DoubleVar()
    Thickness_ = DoubleVar()
    Poisson_ratio = DoubleVar()
    Young_mod = DoubleVar()
    Force = DoubleVar()
    nodes_x = IntVar()
    nodes_y = IntVar()

    Length.set(10)
    Height.set(30)
    Thickness_.set(1)
    Poisson_ratio.set(0.3)
    Young_mod.set(210e9)
    Force.set(10000)
    nodes_x.set(3)
    nodes_y.set(5)

    Entry(root, textvariable=Length).grid(row=1, column=4)
    Entry(root, textvariable=Height).grid(row=2, column=4)
    Entry(root, textvariable=Thickness_).grid(row=3, column=4)
    Entry(root, textvariable=Poisson_ratio).grid(row=4, column=4)
    Entry(root, textvariable=Young_mod).grid(row=5, column=4)
    Entry(root, textvariable=Force).grid(row=6, column=4)
    Entry(root, textvariable=nodes_x).grid(row=7, column=4)
    Entry(root, textvariable=nodes_y).grid(row=8, column=4)
    scvalue1 = DoubleVar()
    scvalue1.set(0)
    scvalue2 = DoubleVar()
    scvalue2.set(0)
    scvalue3 = DoubleVar()
    scvalue3.set(0)
    Entry(root, textvar=scvalue1, font="lucida 8 bold").grid(row=10, column=4, pady=10, padx=10)
    Entry(root, textvar=scvalue2, font="lucida 8 bold").grid(row=11, column=4, pady=10, padx=10)
    Entry(root, textvar=scvalue3, font="lucida 8 bold").grid(row=12, column=4, pady=10, padx=10)

    # disp = ""
    b = Button(root, text="Calculate", command=solve).grid(row=9, column=4, pady=15)
    b = Button(root, text="Close", command=root.destroy).grid(row=13, column=4, pady=15)
    root.mainloop()


if __name__ == "__main__":
    run()
