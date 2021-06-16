from matplotlib import pyplot as plt

with open("ilsrr.dat", "r") as data:
    # Define the max inner steps
    outer_count = 0
    max_inner_steps = 0
    for line in data:
        elem = line.split()
        if elem[0] == "OUTER":
            outer_count += 1
            if outer_count == 1:
                break
        elif elem[0] == "INNER":
            max_inner_steps = int(elem[2])

    outer_steps  = []
    inner_steps  = []
    outer_energy = []
    inner_energy = []
    outer_step   = 0
    for line in data:
        elem = line.split()
        if elem[0] == "OUTER":
            outer_step = int(elem[2])
            outer_steps.append(outer_step * max_inner_steps)
            outer_energy.append(float(elem[5]))
        elif elem[0] == "INNER":
            inner_steps.append(int(elem[2]) + outer_step * max_inner_steps)
            inner_energy.append(float(elem[5]))

print("CALLED")
plt.plot(inner_steps, inner_energy)
plt.show()