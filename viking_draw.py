# Import python libraries.
import matplotlib.pyplot as plt



# Read map data from file.
def calc_xy_from_coords(coords):
    degree_to_km = 37.1133457954
    x = coords[1] * degree_to_km
    y = coords[0] * degree_to_km
    return (x, y)

map_data = []
success_map_data = []

with open("data/map/greenland.dat") as file:
    for line in file:
        data = calc_xy_from_coords((float(line.split()[1]), float(line.split()[0])))
        map_data.append(data)

with open("data/map/greenland_success.dat") as file:
    for line in file:
        data = calc_xy_from_coords((float(line.split()[1]), float(line.split()[0])))
        success_map_data.append(data)



#Do the graphing.
simulation_start = 0
simulation_end = 9

for k in range(1): #2
    if(k == 0):
        simulation_equinox = "sol"
    else:
        simulation_equinox = "equ"

    for l in range(1):
        if(l == 0):
            simulation_crystal = "cal"
        if(l == 1):
            simulation_crystal = "cord"
        if(l == 2):
            simulation_crystal = "tour"

        for m in range(5):
            if(m == 0):
                simulation_h = 0.5
            if(m == 1):
                simulation_h = 1.0
            if(m == 2):
                simulation_h = 2.0
            if(m == 3):
                simulation_h = 3.0
            if(m == 4):
                simulation_h = 6.0

            filename_seed = simulation_equinox + "_" + simulation_crystal + "_" + str(simulation_h)
            print(simulation_equinox + '\t' + simulation_crystal + '\t' + str(simulation_h))

            fig, ax = plt.subplots()
            ax.plot([i[0] for i in map_data], [i[1] for i in map_data], '.', markersize=1, color= 'black')
            ax.plot([i[0] for i in success_map_data], [i[1] for i in success_map_data], color= 'blue')
            ax.set(xlabel='x (km)', ylabel='y (km)', title='Map')
            ax.grid()
            ax.set_xlim(-2550, 750)
            ax.set_ylim(1800, 3150)

            for i in range(simulation_start, simulation_end):

                filename = "simulations/" + filename_seed + '-' + str(i).zfill(5) + ".dat"
                f_in = open(filename, "r")
                col = "red"
                pos_array = []
                with open(filename) as f_in:
                    if(f_in.readline()[1:-1] == "True"):
                        col = "green"
                    for line in f_in:
                        pos_array.append((float(line.split()[4][1:-1]), float(line.split()[5][:-1])))

                ax.plot([i[0] for i in pos_array], [i[1] for i in pos_array], color=col)

            fig.savefig("pics/" + filename_seed + ".png", dpi=300)
            fig.clear()
            plt.close()
