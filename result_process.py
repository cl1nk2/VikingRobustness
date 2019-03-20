#Import python libraries.
import matplotlib.pyplot as plt



result_array = []
with open("result.dat", "r") as f_in:
    for line in f_in:
        result_array.append((line.split()[0], line.split()[1], float(line.split()[2]), bool(line.split()[3] == 'True'), int(line.split()[4]), int(line.split()[5]), float(line.split()[6])))

fig, ax = plt.subplots()

night_prog_dict = {}
not_night_prog_dict = {}
for i in result_array:
    if(i[3]):
        if(i[2] in night_prog_dict):
            night_prog_dict[i[2]] = tuple(map(sum,zip(night_prog_dict[i[2]],(i[6], 1))))
        else:
            night_prog_dict[i[2]] = (i[6], 1)
    else:
        if(i[2] in not_night_prog_dict):
            not_night_prog_dict[i[2]] = tuple(map(sum,zip(not_night_prog_dict[i[2]],(i[6], 1))))
        else:
            not_night_prog_dict[i[2]] = (i[6], 1)
            
for i in night_prog_dict:
    night_prog_dict[i] = night_prog_dict[i][0]/night_prog_dict[i][1]*100
for i in not_night_prog_dict:
    not_night_prog_dict[i] = not_night_prog_dict[i][0]/not_night_prog_dict[i][1]*100

ax.plot([i for i in night_prog_dict], [night_prog_dict[i] for i in night_prog_dict], '-', color="blue")
ax.plot([i for i in not_night_prog_dict], [not_night_prog_dict[i] for i in not_night_prog_dict], '-', color="red")
ax.set(xlabel='Rate of measurement (hours)', ylabel='Success rate (%)', title='Success rate of night progression turned on and off')
ax.grid()
ax.set_xlim(0, 7)
ax.set_ylim(0, 100)
ax.legend(("Night progression turned on", "Night progression turned off"))

fig.savefig("results/night_progression.png", dpi=300)
fig.clear()
plt.close()