from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
first = np.loadtxt(folder / "ca_ow_rdf_half1.xvg", comments=("#", "@"))
second = np.loadtxt(folder / "ca_ow_rdf_half2.xvg", comments=("#", "@"))

assert first.shape == second.shape
assert first.shape[1] == 2
assert np.isfinite(first).all() and np.isfinite(second).all()
assert np.allclose(first[:, 0], second[:, 0])

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

for ax in axes:
    ax.plot(first[:, 0], first[:, 1],
            label="0–999 ps", color="navy")
    ax.plot(second[:, 0], second[:, 1],
            label="1000–2000 ps", color="darkorange", linestyle="--")
    ax.set_xlabel("Ca–O distance (nm)")
    ax.set_ylabel("RDF g(r)")
    ax.grid(alpha=0.2)
    ax.legend()

axes[0].set_xlim(0, 0.6)
axes[0].set_title("Hydration shells")
axes[1].set_xlim(0.22, 0.29)
axes[1].set_title("First peak: enlarged")

fig.suptitle("Ca²⁺ CM — replicate 02 - RDF comparison between trajectory halves")
fig.tight_layout()

output = folder / "ca_ow_rdf_halves.png"
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"图片已保存：{output}")
