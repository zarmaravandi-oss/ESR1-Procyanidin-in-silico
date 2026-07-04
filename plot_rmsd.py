import matplotlib.pyplot as plt

with open("rmsd_complex.xvg") as f:
    lines = [line for line in f if not line.startswith(("#", "@"))]
    time_ps = []
    rmsd = []

    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            time_ps.append(float(parts[0]))
            rmsd.append(float(parts[1]))

time_ns = [t / 1000 for t in time_ps]

plt.figure(figsize=(10, 6))
plt.plot(time_ns, rmsd, color='darkorange', linewidth=2)
plt.xlabel("Time (ns)", fontsize=12)
plt.ylabel("RMSD (nm)", fontsize=12)
plt.title("RMSD of Protein-Ligand Complex", fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig("rmsd_plot_100ns.png", dpi=300)
plt.show()
