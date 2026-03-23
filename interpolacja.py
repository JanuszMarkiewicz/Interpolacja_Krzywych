import matplotlib.pyplot as plt
import numpy as np
import ast
import os
def add_spline_to_plot(x, y, M=150, show_nodes=True):
    """
    Parameters:
    x : list - x-coordinates of nodes
    y : list - y-coordinates of nodes
    M : int - number of points to generate on the curve (default: 150)
    show_nodes : bool - whether to display nodes on the plot (default: True)
    """

    def compute_second_derivatives(t, y): # from NIFS3 algorithm
        """
        Function to compute the second derivatives M_k (s'') for cubic splines
        according to the algorithm.

        Parameters:
        t : list - interpolation nodes (x_k)
        y : list - function values at the nodes (f(x_k))

        Returns:
        M : list - second derivatives at the nodes (M_k)
        """
        n = len(t)
        h = [t[i] - t[i-1] for i in range(1, len(t))]
        lambd = [h[i] / (h[i] + h[i+1]) for i in range(n-2)] 

        d = [6 * (((y[i+1] - y[i]) / h[i]) - ((y[i] - y[i-1]) / h[i-1])) / (t[i+1] - t[i-1]) for i in range(1, n-1)]

        p = [0] * n
        q = [0] * n
        u = [0] * n

        q[0] = u[0] = 0
        for k in range(1, n-2):
            p[k] = lambd[k-1] * q[k-1] + 2
            q[k] = (lambd[k-1] - 1) / p[k]
            u[k] = (d[k-1] - lambd[k-1] * u[k-1]) / p[k]

        M = [0] * (n + 1)
        M[n-1] = u[n-1]
        for k in range(n-3, 0, -1):
            M[k] = u[k] + q[k] * M[k+1]

        return M, h

    def evaluate_spline(t, y, M, h, u): # from NIFS3 theorem
        """
        Function to evaluate spline.
        
        Parameters:
        t : list - parameter values of the nodes
        y : list - values of the function at the nodes
        M : list - second derivatives at the nodes
        h : list - intervals between nodes
        u : float - parameter value at which to evaluate the spline

        Returns:
        float - interpolated value of the spline at u
        """
        n = len(t)
        for i in range(n-1):
            if t[i] <= u <= t[i+1]:
                dx = u - t[i]
                term1 = (M[i] * (t[i+1] - u)**3) / (6 * h[i])
                term2 = (M[i+1] * (dx**3)) / (6 * h[i])
                term3 = ((y[i] - (M[i] * h[i]**2) / 6) * (t[i+1] - u)) / h[i]
                term4 = ((y[i+1] - M[i+1] * h[i]**2 / 6) * dx) / h[i]
                return term1 + term2 + term3 + term4
        return None


    N = len(x)
    t = [i/(N-1) for i in range(N)]  # t=[k/(N-1)]

    M_x, h_x = compute_second_derivatives(t, x)
    M_y, h_y = compute_second_derivatives(t, y)

    u_values = np.linspace(0, 1, M)  # u=[i/M], i=0,1,...,M-1

    spline_x = [evaluate_spline(t, x, M_x, h_x, u) for u in u_values]
    spline_y = [evaluate_spline(t, y, M_y, h_y, u) for u in u_values]

    # interchangable
    plt.plot(spline_x, spline_y, color = "red") # constant red color
    #plt.plot(spline_x, spline_y) # different curves separated by color

    if show_nodes:
        plt.scatter(x, y)

def show_plot():
    """
    Function that draws all added curves and nodes
    """
    plt.title("NIFS3")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.show()

def get_from_file(file_name):
    """
    Reads x and y data from a file in the 'data' folder and returns them as lists.
    """
    file_path = os.path.join("data", file_name)
    
    with open(file_path, "r") as file:
        data = file.read()
    
    x_str = data.split("x:=")[1].split(",\n")[0].strip()
    y_str = data.split("y:=")[1].split(",\n")[0].strip()

    x = ast.literal_eval(x_str)
    y = ast.literal_eval(y_str)

    return x, y

def save_to_file(data, filename, M):
    """
    Converts (x, y) data into x, y, t, and u arrays, then appends them to a file in the specified format.

    Parameters:
        data (tuple): Tuple containing x and y coordinate lists.
        filename (str): Name of the file to save the data.
        M (int): Number of subdivisions for u values.
    """
    x, y = data

    N = len(x)
    t = [i / (N - 1) for i in range(N)]
    u_values = np.linspace(0, 1, M).tolist()

    with open(filename, "a") as file:
        file.write(f"x:={x}, ")
        file.write(f"y:={y}, ")
        file.write(f"t:={t}, ")
        file.write(f"u:={u_values}\n\n")

def save_plot(filename="plot.jpg", dpi=300, linewidth_scale=0.6):
    """
    Saves the current plot to a file with high resolution and consistent line thickness.

    Parameters:
    filename : str - Name of the file to save the plot (default: "plot.png")
    dpi : int - Dots per inch (resolution) for the saved plot 
    linewidth_scale : float - Scale factor for line width 
    """
    ax = plt.gca()
    
    for line in ax.get_lines():
        line.set_linewidth(line.get_linewidth() * linewidth_scale)
    
    ax.set_aspect('equal', adjustable='datalim')
    plt.gcf().set_size_inches(16, 8)
    
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    print(f"Plot saved as {filename} with resolution {dpi} dpi.")

# --------------------------------------------------------------------------------------------------------- #
#                              Here is the description of how each NIFS3 is specified
# --------------------------------------------------------------------------------------------------------- #

boolean = False
filename = "konkurs-I-351009-dane.txt"

# each is a different NIFS3:

# P
x1, y1 = get_from_file("P.txt") # get pixel coordinates from file(each curve has a designated data file), in this case P.txt holds pixel values for drawing 'P'
add_spline_to_plot(x1, y1, M=60, show_nodes=boolean) # calculate
save_to_file((x1,y1), filename, M=60) # save data to file

# CONTINUE ACCORDINGLY FOR EACH CURVE:

# W
x2, y2 = get_from_file("W.txt")
add_spline_to_plot(x2, y2, M=70, show_nodes=boolean)
save_to_file((x2, y2), filename, M=70)

# O
x3, y3 = get_from_file("O.txt")
add_spline_to_plot(x3, y3, M=50, show_nodes=boolean)
save_to_file((x3, y3), filename, M=50)

# +
x4, y4 = get_from_file("plus11.txt")
add_spline_to_plot(x4, y4, M=5, show_nodes=boolean)
save_to_file((x4, y4), filename, M=5)

x5, y5 = get_from_file("plus12.txt")
add_spline_to_plot(x5, y5, M=5, show_nodes=boolean)
save_to_file((x5, y5), filename, M=5)

# +
x6, y6 = get_from_file("plus21.txt")
add_spline_to_plot(x6, y6, M=25, show_nodes=boolean)
save_to_file((x6, y6), filename, M=25)

# t1, t2, t3
x7, y7 = get_from_file("t1.txt")
add_spline_to_plot(x7, y7, M=5, show_nodes=boolean)
save_to_file((x7, y7), filename, M=5)

x8, y8 = get_from_file("t2.txt")
add_spline_to_plot(x8, y8, M=2, show_nodes=boolean)
save_to_file((x8, y8), filename, M=2)

x9, y9 = get_from_file("t3.txt")
add_spline_to_plot(x9, y9, M=2, show_nodes=boolean)
save_to_file((x9, y9), filename, M=2)

# o
x10, y10 = get_from_file("o.txt")
add_spline_to_plot(x10, y10, M=30, show_nodes=boolean)
save_to_file((x10, y10), filename, M=30)

# n
x11, y11 = get_from_file("n.txt")
add_spline_to_plot(x11, y11, M=35, show_nodes=boolean)
save_to_file((x11, y11), filename, M=35)

# aj, aj2
x12, y12 = get_from_file("aj.txt")
add_spline_to_plot(x12, y12, M=60, show_nodes=boolean)
save_to_file((x12, y12), filename, M=60)

x13, y13 = get_from_file("aj2.txt")
add_spline_to_plot(x13, y13, M=3, show_nodes=boolean)
save_to_file((x13, y13), filename, M=3)

# le
x14, y14 = get_from_file("le.txt")
add_spline_to_plot(x14, y14, M=50, show_nodes=boolean)
save_to_file((x14, y14), filename, M=50)

# p
x15, y15 = get_from_file("p.txt")
add_spline_to_plot(x15, y15, M=50, show_nodes=boolean)
save_to_file((x15, y15), filename, M=50)

# rzy
x16, y16 = get_from_file("rzy.txt")
add_spline_to_plot(x16, y16, M=70, show_nodes=boolean)
save_to_file((x16, y16), filename, M=70)

# j21, j22
x17, y17 = get_from_file("j21.txt")
add_spline_to_plot(x17, y17, M=8, show_nodes=boolean)
save_to_file((x17, y17), filename, M=8)

x18, y18 = get_from_file("j22.txt")
add_spline_to_plot(x18, y18, M=3, show_nodes=boolean)
save_to_file((x18, y18), filename, M=3)

# e
x19, y19 = get_from_file("e.txt")
add_spline_to_plot(x19, y19, M=50, show_nodes=boolean)
save_to_file((x19, y19), filename, M=50)

# zy
x20, y20 = get_from_file("zy.txt")
add_spline_to_plot(x20, y20, M=60, show_nodes=boolean)
save_to_file((x20, y20), filename, M=60)

# k
x21, y21 = get_from_file("k.txt")
add_spline_to_plot(x21, y21, M=35, show_nodes=boolean)
save_to_file((x21, y21), filename, M=35)

# pr
x22, y22 = get_from_file("pr.txt")
add_spline_to_plot(x22, y22, M=80, show_nodes=boolean)
save_to_file((x22, y22), filename, M=80)

# o2
x23, y23 = get_from_file("o2.txt")
add_spline_to_plot(x23, y23, M=45, show_nodes=boolean)
save_to_file((x23, y23), filename, M=45)

# gr
x24, y24 = get_from_file("gr.txt")
add_spline_to_plot(x24, y24, M=65, show_nodes=boolean)
save_to_file((x24, y24), filename, M=65)

# a
x25, y25 = get_from_file("a.txt")
add_spline_to_plot(x25, y25, M=35, show_nodes=boolean)
save_to_file((x25, y25), filename, M=35)

# m
x26, y26 = get_from_file("m.txt")
add_spline_to_plot(x26, y26, M=40, show_nodes=boolean)
save_to_file((x26, y26), filename, M=40)

# o3
x27, y27 = get_from_file("o3.txt")
add_spline_to_plot(x27, y27, M=45, show_nodes=boolean)
save_to_file((x27, y27), filename, M=45)

# w
x28, y28 = get_from_file("w.txt")
add_spline_to_plot(x28, y28, M=50, show_nodes=boolean)
save_to_file((x28, y28), filename, M=50)

# an
x29, y29 = get_from_file("an.txt")
add_spline_to_plot(x29, y29, M=90, show_nodes=boolean)
save_to_file((x29, y29), filename, M=90)

# i1, i2
x30, y30 = get_from_file("i1.txt")
add_spline_to_plot(x30, y30, M=15, show_nodes=boolean)
save_to_file((x30, y30), filename, M=15)

x31, y31 = get_from_file("i2.txt")
add_spline_to_plot(x31, y31, M=3, show_nodes=boolean)
save_to_file((x31, y31), filename, M=3)

# a2
x32, y32 = get_from_file("a2.txt")
add_spline_to_plot(x32, y32, M=40, show_nodes=boolean)
save_to_file((x32, y32), filename, M=40)

# oko1, oko2
x33, y33 = get_from_file("oko1.txt")
add_spline_to_plot(x33, y33, M=20, show_nodes=boolean)
save_to_file((x33, y33), filename, M=20)

x34, y34 = get_from_file("oko2.txt")
add_spline_to_plot(x34, y34, M=13, show_nodes=boolean)
save_to_file((x34, y34), filename, M=13)

# nos
x35, y35 = get_from_file("nos.txt")
add_spline_to_plot(x35, y35, M=12, show_nodes=boolean)
save_to_file((x35, y35), filename, M=12)

# usta
x36, y36 = get_from_file("usta.txt")
add_spline_to_plot(x36, y36, M=10, show_nodes=boolean)
save_to_file((x36, y36), filename, M=10)

save_plot("konkurs-I-351009.jpg", dpi=800, linewidth_scale=0.6) # save the plot
#show_plot()


# --------------------------------------------------------------------------------------------------------- #

def process_file(file_path, output_file): # summary file generator
    num_records = 0
    num_interpolation_points_total = 0
    sum_u_lengths_total = 0
    
    with open(file_path, 'r') as f:
        data = f.read().strip()

    records = data.split("\n\n")  
    
    for record in records:
        num_records += 1
        
        t_data = record.split("t:=")[1].split(", u:=")[0].strip()
        u_data = record.split("u:=")[1].strip()
            
        t = ast.literal_eval(t_data)
        u = ast.literal_eval(u_data)

        
        num_interpolation_points_total += len(t)
        sum_u_lengths_total += len(u)

    
    with open(output_file, 'w') as out_file:
        out_file.write(f"{num_records}, {num_interpolation_points_total}, {sum_u_lengths_total}\n")

process_file("konkurs-I-351009-dane.txt","konkurs-I-351009-podsumowanie.txt") # sum up values from input file and save in output file
