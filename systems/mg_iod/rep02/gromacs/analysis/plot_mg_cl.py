from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
data = np.loadtxt(folder / "mg_cl_distances.xvg", comments=("#", "@"))

assert data.ndim == 2 and data.shape[1] == 3
assert np.isfinite(data).all()
assert np.all(np.diff(data[:, 0]) > 0)
assert np.all(data[:, 1:] >= 0)

time_ns = data[:, 0] / 1000

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, sharey=True)

for column, ax, color in zip((1, 2), axes, ("navy", "darkorange")):
    distance = data[:, column]
    closest = np.argmin(distance)

    ax.plot(time_ns, distance, color=color, linewidth=0.9)
    ax.plot(time_ns[closest], distance[closest], "o", color="black",
            label=f"Minimum: {distance[closest]:.3f} nm")
    ax.set_ylabel(f"Mg–Cl{column} distance (nm)")
    ax.set_ylim(0, data[:, 1:].max() * 1.05)
    ax.grid(alpha=0.2)
    ax.legend()

axes[-1].set_xlabel("Time (ns)")
fig.suptitle("Mg²⁺ IOD — replicate 02 - distances to the two chloride ions")
fig.tight_layout()

output = folder / "mg_cl_distances.png"
assert not output.exists(), output
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"图片已保存：{output}")
