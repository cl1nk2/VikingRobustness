#Import python libraries.
import matplotlib.pyplot as plt



result_array = []
with open("result.dat", "r") as f_in:
    for line in f_in:
        result_array.append((line.split()[0], line.split()[1], float(line.split()[2]), bool(line.split()[3]), int(line.split()[4]), int(line.split()[5]), float(line.split()[6])))

fig, ax = plt.subplots()

for i in result_array:
    if(i[0] == "equ"):
        if(i[1] == "cord"):
            if(not(i[3])):
                if(i[4] == 2):
                    if(i[5] == 1):
                        ax.plot(i[2], i[6], '.', markersize=1, color= 'black')

ax.set(xlabel='parameter1', ylabel='Success', title='Success rate of parameter1')
ax.grid()
#ax.set_xlim(0, 1)
#ax.set_ylim(0, 1)

fig.savefig("results/parameter1.png", dpi=300)
fig.clear()
plt.close()