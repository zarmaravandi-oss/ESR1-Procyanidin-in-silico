import matplotlib.pyplot as plt

with open("rg_complex.xvg") as f:
    lines = [line for line in f if not line.startswith(("#", "@"))]
    time_ps = []
    rg = []

    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            time_ps.append(float(parts[0]))
            rg.append(float(parts[1]))

time_ns = [t / 1000 for t in time_ps]

plt.plot(time_ns, rg)
plt.xlabel("Time (ns)")
plt.ylabel("Radius of Gyration (nm)")
plt.title("RG of Protein-Ligand Complex")
plt.grid(True)
plt.show()
