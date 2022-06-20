import os
import copy
import math

import matplotlib.pyplot

import pyopal.objects.parser
import PyOpal.elements.vertical_ffa_magnet

print("Running vertical_ffa_magnet test")

def get_div_b(magnet, pos):
    dx = 0.1
    db = [0., 0., 0.]
    for i in range(3):
        my_pos = copy.deepcopy(pos)
        my_pos[i] += dx
        db_plus = magnet.get_field_value(my_pos[0], my_pos[1], my_pos[2], my_pos[3])[i+1]
        my_pos[i] -= 2*dx
        db_minus = magnet.get_field_value(my_pos[0], my_pos[1], my_pos[2], my_pos[3])[i+1]
        db[i] = (db_plus-db_minus)/2/dx
    return sum(db)/10*1e-3 # kGauss/mm to T/m

def get_curl_1d(magnet, pos, a, b):
    # db_b/dx_a - db_a/dx_b
    dx = 0.1
    db = []
    for i, j in (a, b), (b, a):
        my_pos = copy.deepcopy(pos)
        my_pos[i] += dx
        db_plus = magnet.get_field_value(my_pos[0], my_pos[1], my_pos[2], my_pos[3])[j+1]
        my_pos[i] -= 2*dx
        db_minus = magnet.get_field_value(my_pos[0], my_pos[1], my_pos[2], my_pos[3])[j+1]
        db.append((db_plus-db_minus)/2/dx)
    return db[0] - db[1]


def get_curl_b_mag(magnet, pos):
    curl_mag = get_curl_1d(magnet, pos, 3, 2)**2+\
               get_curl_1d(magnet, pos, 3, 1)**2+\
               get_curl_1d(magnet, pos, 1, 2)**2
    curl_mag = curl_mag**0.5
    print (curl_mag)
    return curl_mag

def draw_polynomial(axes, cell_length, n_cells, centre_x, centre_y):
    x_list, y_list = [centre_x], [-cell_length/2.0+centre_y]
    theta = 0
    for i in range(n_cells):
        x_list.append(x_list[-1]-cell_length*math.sin(theta))
        y_list.append(y_list[-1]+cell_length*math.cos(theta))
        theta += math.pi*2.0/n_cells
    axes.plot(x_list, y_list)

def plot(magnet, x_range, y_range, max_b, max_db):
    x_list, y_list, z1_list, z2_list, z3_list, z4_list = [], [], [], [], [], []
    n_bins = 201
    for i in range(n_bins):
        for j in range(n_bins):
            x = i*(x_range[1]-x_range[0])/(n_bins-1)+x_range[0]
            z = j*(y_range[1]-y_range[0])/(n_bins-1)+y_range[0]
            field = list(magnet.get_field_value(x, 0.0, z, 0.0))
            for k in range(1, 4):
                field[k] = field[k]/10.0 # kGauss to T
            b = (field[1]**2+field[2]**2+field[3]**2)**0.5
            if b > max_b:
                b = 0
            by = abs(field[2])
            if by > max_b:
                by = 0
            div_b = abs(get_div_b(magnet, [x, 0.0, z, 0.0]))
            curl_b = get_curl_b_mag(magnet, [x, 0.0, z, 0.0])
            if abs(div_b) > max_db:
                div_b = 0 
            x_list.append(x)
            y_list.append(z)
            z1_list.append(b)
            z2_list.append(by)
            z3_list.append(div_b)
            z4_list.append(curl_b)
    hrange = [[min(x_list)-1.0, max(x_list)+1.0], [min(y_list)-5.0, max(y_list)+5.0]]
    figure = matplotlib.pyplot.figure(figsize=(20, 10))

    axes = figure.add_subplot(2, 2, 1)
    hist2d = axes.hist2d(x_list, y_list, n_bins, hrange, weights=z1_list)
    axes.get_figure().colorbar(hist2d[3], ax=axes)
    axes.set_title("B$_{tot}$ [T]")
    draw_polynomial(axes, 2800, 10, 0.0, sum(y_range)/2.0+530.0)

    axes = figure.add_subplot(2, 2, 2)
    hist2d = axes.hist2d(x_list, y_list, n_bins, hrange, weights=z2_list)
    axes.get_figure().colorbar(hist2d[3], ax=axes)
    axes.set_title("B$_y$ [T]")
    draw_polynomial(axes, 2800, 10, 0.0, sum(y_range)/2.0+530.0)

    axes = figure.add_subplot(2, 2, 3)
    hist2d = axes.hist2d(x_list, y_list, n_bins, hrange, weights=z3_list, norm=matplotlib.colors.LogNorm())
    axes.get_figure().colorbar(hist2d[3], ax=axes)
    axes.set_title("|div(b)| [T/m]")
    draw_polynomial(axes, 2800, 10, 0.0, sum(y_range)/2.0+530.0)

    axes = figure.add_subplot(2, 2, 4)
    hist2d = axes.hist2d(x_list, y_list, n_bins, hrange, weights=z4_list, norm=matplotlib.colors.LogNorm())
    axes.get_figure().colorbar(hist2d[3], ax=axes)
    axes.set_title("|curl(b)| [T/m]")
    draw_polynomial(axes, 2800, 10, 0.0, sum(y_range)/2.0+530.0)
    figure.savefig("field_2d_"+magnet.end_field_model+".png")

    matplotlib.pyplot.show(block=False)

def plot_maxwell(magnet, pos, order_list):
    div_b_list = []    
    curl_b_list = []    
    for i, order in enumerate(order_list):
        magnet.max_horizontal_power = order
        div_b_list.append(abs(get_div_b(magnet, pos)))
        curl_b_list.append(get_curl_b_mag(magnet, pos))

    figure = matplotlib.pyplot.figure(figsize=(20, 10))
    axes = figure.add_subplot(1, 1, 1)
    axes.set_title(magnet.end_field_model)
    axes.plot(order_list, div_b_list, label="|Div(B)| [T/m]")
    axes.plot(order_list, curl_b_list, label="|Curl(B)| [T/m]")
    axes.set_xlabel("Order")
    axes.set_ylabel("Derivative [T/m]")
    axes.semilogy()
    axes.legend()
    figure.savefig("maxwell_"+magnet.end_field_model+".png")

def main():
    magnet = PyOpal.vertical_ffa_magnet.VerticalFFAMagnet()
    magnet.b0 = 1.0
    magnet.field_index = 1.31
    magnet.max_horizontal_power = 4
    magnet.centre_length = 0.5
    magnet.end_length = 0.15
    magnet.end_field_model = "arctan"
    magnet.width = 2*math.sin(math.radians(36))*2.8
    magnet.height_neg_extent = 1.0
    magnet.height_pos_extent = 1.0
    magnet.bb_length = 12.0
    magnet.enge_parameters = [1.0, 2.0, 3.0]
    magnet.bb_start_position = [0.0, 0.0, 0.0]
    magnet.bb_start_normal = [0.0, 0.0, 1.0]
    magnet.bb_end_position = [0.0, 0.0, 12.0]
    magnet.bb_end_normal = [0.5, 0.0, 0.5]

    print("VFFA lookup")
    magnet.get_field_value(100.0, 0.0, 6000.0-250.0, 0.0)
    plot(magnet, [-4000.0, 3000.0], [500.0, 12000.0], 5, 0.1)
    #plot_maxwell(magnet, [100.0, 100.0, 4000.0, 0.0], [i for i in range(11)])

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Done")
