output = open("greenland_success.dat", "w+")

with open("map.dat") as file:
    for line in file:
        longitude = float(line.split()[1])
        latitude = float(line.split()[0])
        if(latitude > -45.9 and latitude < -39.5 and longitude < 65.5 and longitude > 59.5):
            output.write(line)