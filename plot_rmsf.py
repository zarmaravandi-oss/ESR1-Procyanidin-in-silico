import matplotlib.pyplot as plt

with open("rmsf.xvg") as f:
    lines = [line for line in f if not line.startswith(("#", "@"))]
    residue = []
    fluctuation = []

    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            residue.append(int(float(parts[0])))
            fluctuation.append(float(parts[1]))

plt.figure(figsize=(10, 6))
plt.plot(residue, fluctuation, color='mediumslateblue', linewidth=2)
plt.xlabel("Residue Number", fontsize=12)
plt.ylabel("RMSF (nm)", fontsize=12)
plt.title("RMSF per Residue in Protein-Ligand Complex", fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig("rmsf_plot.png", dpi=300)
plt.show()
