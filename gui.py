from tkinter import *

root = Tk()
root.geometry("644x344")


def getvals():
    print(length.get())
    print(height.get())
    print(thickness.get())


Label(root, text="Welcome to FEM World", font="comicsansms 13 bold", pady=15).grid(row=0, column=4)
# Text for our form
Label(root, text = "Length").grid(row = 1, column = 2)
Label(root, text = "Height").grid(row = 2, column = 2)
Label(root, text = "Thickness").grid(row = 3, column = 2)
Label(root, text = "Poisson Ratio").grid(row = 4, column = 2)
Label(root, text = "Young Modulus").grid(row = 5, column = 2)
Label(root, text = "No of nodes in x direction").grid(row = 6, column = 2)
Label(root, text = "No of nodes in y direction").grid(row = 7, column = 2)

length = DoubleVar()
height = DoubleVar()
thickness = DoubleVar()
poisson_ratio = DoubleVar()
E = DoubleVar()
nodes_x = IntVar()
nodes_y = IntVar()

Entry(root, textvariable = length).grid(row = 1, column = 4)
Entry(root, textvariable = height).grid(row = 2, column = 4)
Entry(root, textvariable = thickness).grid(row = 3, column = 4)
Entry(root, textvariable = poisson_ratio).grid(row = 4, column = 4)
Entry(root, textvariable = E).grid(row = 5, column = 4)
Entry(root, textvariable = nodes_x).grid(row = 6, column = 4)
Entry(root, textvariable = nodes_y).grid(row = 7, column = 4)

Button(text = "Submit", command = getvals).grid(row = 8, column = 4, pady = 15)
root.mainloop()