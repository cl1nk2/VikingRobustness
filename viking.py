# Import python libraries.
import random as rd
import numpy as np



# Read elevation data from files.
equ_elevation_table = {}
sol_elevation_table = {}

with open("data/elevation/elevation_Bergen_equ.dat") as file:
    next(file)
    for line in file:
        equ_elevation_table[float(line.split()[0])] = float(line.split()[1])

with open("data/elevation/elevation_Bergen_sol.dat") as file:
    next(file)
    for line in file:
        sol_elevation_table[float(line.split()[0])] = float(line.split()[1])



# Read error data from files.
cal_equ_am_table = {}
cal_equ_pm_table = {}
cal_sol_am_table = {}
cal_sol_pm_table = {}

with open("data/n_error/cal_equ_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cal_equ_am_table:
            cal_equ_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cal_equ_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cal_equ_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cal_equ_pm_table:
            cal_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cal_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cal_sol_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cal_sol_am_table:
            cal_sol_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cal_sol_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cal_sol_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cal_sol_pm_table:
            cal_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cal_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

cord_equ_am_table = {}
cord_equ_pm_table = {}
cord_sol_am_table = {}
cord_sol_pm_table = {}

with open("data/n_error/cord_equ_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cord_equ_am_table:
            cord_equ_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cord_equ_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cord_equ_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cord_equ_pm_table:
            cord_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cord_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cord_sol_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cord_sol_am_table:
            cord_sol_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cord_sol_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/cord_sol_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in cord_sol_pm_table:
            cord_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            cord_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

tour_equ_am_table = {}
tour_equ_pm_table = {}
tour_sol_am_table = {}
tour_sol_pm_table = {}

with open("data/n_error/tour_equ_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in tour_equ_am_table:
            tour_equ_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            tour_equ_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/tour_equ_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in tour_equ_pm_table:
            tour_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            tour_equ_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/tour_sol_am.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in tour_sol_am_table:
            tour_sol_am_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            tour_sol_am_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

with open("data/n_error/tour_sol_pm.csv") as file:
    next(file)
    for line in file:
        if (int(line.split()[0]), int(line.split()[1])) in tour_sol_pm_table:
            tour_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))].append(int(line.split()[2]))
        else:
            tour_sol_pm_table[(int(line.split()[0]), int(line.split()[1]))] = [int(line.split()[2])]

error_tables = {("cal", "equ", "am") : cal_equ_am_table,
                ("cal", "equ", "pm") : cal_equ_pm_table,
                ("cal", "sol", "am") : cal_sol_am_table,
                ("cal", "sol", "pm") : cal_sol_pm_table,
                ("cord", "equ", "am") : cord_equ_am_table,
                ("cord", "equ", "pm") : cord_equ_pm_table,
                ("cord", "sol", "am") : cord_sol_am_table,
                ("cord", "sol", "pm") : cord_sol_pm_table,
                ("tour", "equ", "am") : tour_equ_am_table,
                ("tour", "equ", "pm") : tour_equ_pm_table,
                ("tour", "sol", "am") : tour_sol_am_table,
                ("tour", "sol", "pm") : tour_sol_pm_table}



# Read map data from file.
def calc_xy_from_coords(coords):
    degree_to_km = 37.1133457954
    x = coords[1] * degree_to_km
    y = coords[0] * degree_to_km
    return (x, y)

success_map_data = []

with open("data/map/greenland_success.dat") as file:
    for line in file:
        data = calc_xy_from_coords((float(line.split()[1]), float(line.split()[0])))
        success_map_data.append(data)



# Define function for calculating.
def calc_meas_time(h):
    return int(round(rd.uniform(h-h/6.0, h+h/6.0)*60))

def calc_wake_up_time(equ):
    if(equ):
        return rd.uniform(5.0, 6.0)
    else:
        return rd.uniform(2.0, 3.0)

def calc_elevation(equ, hour):
	elevation = 0.0
	hour = int(hour)

	if(equ):
		if(hour > 5 and hour < 18):
			elevation = equ_elevation_table[hour]
		else:
			elevation = 0.0

	else:
		if(hour > 2 and hour < 21):
			elevation = sol_elevation_table[hour]
		else:
			elevation = 0.0

	return elevation

def calc_cloudiness(cloudiness, cloud_med, cloud_dev):
	delta = rd.gauss(cloud_med, cloud_dev)
	cloudiness = int(round(cloudiness + delta))
	if(cloudiness > 8):
		cloudiness = 8
	if(cloudiness < 0):
		cloudiness = 0
	return cloudiness

def calc_n_error(equ, cal, cord, tour, hour, elevation_double, cloudiness):
    elevation = int(np.ceil(elevation_double / 5.0) * 5.0)
    if(elevation < 5):
        elevation = 5
    if(equ):
        if(elevation >= 25):
            elevation = 25
    else:
        if(elevation >= 50):
            elevation = 50

    n_error = 0

    str1 = ""
    str2 = ""
    str3 = ""

    if(cal):
        str1 = "cal"
    if(cord):
        str1 = "cord"
    if(tour):
        str1 = "tour"

    if(equ):
        str2 = "equ"
    else:
        str2 = "sol"

    if(hour < 12.0):
        str3 = "am"
    else:
        str3 = "pm"

    weather_situation = int(round(rd.uniform(0, len(error_tables[(str1, str2, str3)][(elevation, cloudiness)]) - 1)))

    n_error = error_tables[(str1, str2, str3)][(elevation, cloudiness)][weather_situation]

    return n_error

def calc_velocity(vel_avg, vel_std, vel_max, vel, n_error):
    vel_abs = rd.gauss(vel_avg, vel_std)
    while(vel_abs > vel_max):
        vel_abs = rd.gauss(vel_avg, vel_std)

    phi = np.pi - np.pi / 180.0 * float(n_error)

    vel = (np.cos(phi) * vel_abs, np.sin(phi) * vel_abs)

    return vel



#Hit-box method for finding successful route.
def is_it_a_hit(pos_current, min_distance):
    for pos_target in success_map_data:
        delta_x = pos_target[0] - pos_current[0]
        if(delta_x <= min_distance):
            delta_y = pos_target[1] - pos_current[1]
            if(delta_y <= min_distance):
                distance = np.sqrt(np.power(delta_x, 2) + np.power(delta_y, 2))
                if(distance <= min_distance):
                    return True

    return False

def calc_min_distance():
    return rd.uniform(0.0, 128.23)



#Step, and route functions.
def step(equ, cal, cord, tour, h, vel_avg, vel_std, vel_max, cloud_med, cloud_dev, night_navigation,
            t, meas_time, time_until_meas, wake_up_time, change_wake_up_time,
            elevation, cloudiness, n_error, pos, vel, hit, min_distance):
    pos = (pos[0] + vel[0] / 60.0, pos[1] + vel[1] / 60.0)

    if(t % 60 == 0):
        cloudiness = calc_cloudiness(cloudiness, cloud_med, cloud_dev)

    hour = (t / 60) % 24
    if(equ):
        if(hour > wake_up_time and hour < 18.0):
            if(change_wake_up_time == False):
                change_wake_up_time = True

            if(t > (138000 / vel_max)):
                if(t % 60 == 0):
                    min_distance = calc_min_distance()
                hit = is_it_a_hit(pos, min_distance)

            if(time_until_meas == 0):
                elevation = calc_elevation(equ, int(round(hour)))

                n_error = 0
                n_error = calc_n_error(equ, cal, cord, tour, hour, elevation, cloudiness)
                vel = calc_velocity(vel_avg, vel_std, vel_max, vel, n_error)

                meas_time = calc_meas_time(h)
                time_until_meas = meas_time

        else:
            elevation = 0.0
            n_error = 0
            if(not(night_navigation)):
                vel = (0.0, 0.0)
            if(change_wake_up_time == True):
                wake_up_time = calc_wake_up_time(equ)
                change_wake_up_time = False

    else:
        if(hour > wake_up_time and hour < 21.0):
            if(change_wake_up_time == False):
                change_wake_up_time = True

            if(t > (138000 / vel_max)):
		if(t % 60 == 0):
                    min_distance = calc_min_distance()
                hit = is_it_a_hit(pos, min_distance)

            if(time_until_meas == 0):
                elevation = calc_elevation(equ, int(round(hour)))

                n_error = 0
                n_error = calc_n_error(equ, cal, cord, tour, hour, elevation, cloudiness)
                vel = calc_velocity(vel_avg, vel_std, vel_max, vel, n_error)

                meas_time = calc_meas_time(h)
                time_until_meas = meas_time

        else:
            elevation = 0.0
            n_error = 0
            if(not(night_navigation)):
                vel = (0.0, 0.0)
            if(change_wake_up_time == True):
                wake_up_time = calc_wake_up_time(equ)
                change_wake_up_time = False

    t += 1
    if(time_until_meas > 0):
        time_until_meas -= 1

    return t, meas_time, time_until_meas, wake_up_time, change_wake_up_time, elevation, cloudiness, n_error, pos, vel, hit, min_distance

def route(equ_or_sol, crystal, days, h,
            vel_avg=11.0, vel_std=2.0, vel_max=21.0, cloud_med=0, cloud_dev=2,
            night_navigation=False):
    equ = False
    if(equ_or_sol == "equ"):
        equ = True

    cal = False
    cord = False
    tour = False
    if(crystal == "cal"):
        cal = True
    if(crystal == "cord"):
        cord = True
    if(crystal == "tour"):
        tour = True

    pos = calc_xy_from_coords((61.0, 5.3))
    vel = (-vel_avg, 0.0)

    min_distance = 0.0

    elevation = 0.0
    cloudiness = int(round(rd.uniform(0, 8)))
    n_error = 0

    t = 0
    meas_time = 0
    time_until_meas = 0
    wake_up_time = calc_wake_up_time(equ)
    change_wake_up_time = False

    cloudiness_arr = []
    hit = False

    while (t < (24 * 60 * days)) and not(hit):
    	cloudiness_arr.append(cloudiness)
    	t, meas_time, time_until_meas, wake_up_time, change_wake_up_time, elevation, cloudiness, n_error, pos, vel, hit, min_distance = step(equ, cal, cord, tour, h, vel_avg, vel_std, vel_max, cloud_med, cloud_dev, night_navigation, t, meas_time, time_until_meas, wake_up_time, change_wake_up_time, elevation, cloudiness, n_error, pos, vel, hit, min_distance)
    cloudiness_arr.append(cloudiness)
    
    okta_avg = np.average(cloudiness_arr)
    okta_stdev = np.std(cloudiness_arr)
    
    #TODO: Do the okta_avg, okta_stdev calculations.

    return okta_avg, okta_stdev, hit



#Do the calculations.

print("simulation_equinox,simulation_crystal,simulation_h,cloud_med,cloud_dev,night_navigation,okta_avg,okta_stdev,simul_no,success")
#f_out = open("result.dat", "w+")
simul_no = 0
for i in range(1000000):
    
    simulation_equinox = rd.choice(['equ', 'sol'])
    simulation_crystal = rd.choice(['cal', 'cord', 'tour'])
    if(simulation_equinox == 'equ'):
        days = 16
    if(simulation_equinox == 'sol'):
        days = 13
    simulation_h = rd.uniform(0.5, 6.0)
    cloud_med = rd.uniform(-1, 1)
    cloud_dev = rd.uniform(1, 4)
    night_navigation = rd.choice([False, True])

    #print(simulation_equinox + '\t' + simulation_crystal + '\t' + str(simulation_h) + '\t' + str(cloud_med) + '\t' + str(cloud_dev) + '\t' + str(int(night_navigation)) + '\t' + str(simul_no))
    okta_avg, okta_stdev, hit = route(simulation_equinox, simulation_crystal, days, simulation_h, cloud_med=cloud_med, cloud_dev=cloud_dev, night_navigation=night_navigation)
    print(simulation_equinox + ',' + simulation_crystal + ',' + str(simulation_h) + ',' + str(cloud_med) + ',' + str(cloud_dev) + ',' + str(night_navigation) + ',' + str(okta_avg) + ',' + str(okta_stdev) + ',' + str(simul_no) + ',' + str(int(hit)))
    #f_out.write(simulation_equinox + '\t' + simulation_crystal + '\t' + str(simulation_h) + '\t' + str(cloud_med) + '\t' + str(cloud_dev) + '\t' + str(int(night_navigation)) + '\t' + str(okta_avg) + '\t' + str(okta_stdev) + '\t' + str(simul_no) + '\t' + str(int(hit)) + '\n')
    
    simul_no = simul_no + 1
#f_out.close()		 
