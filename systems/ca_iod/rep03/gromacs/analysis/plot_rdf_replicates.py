from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
root = folder.parents[2]
paths = [
    root / "gromacs" / "analysis" / "ca_ow_rdf.xvg",
    root / "rep02" / "gromacs" / "analysis" / "ca_ow_rdf.xvg",
    folder / "ca_ow_rdf.xvg",
]
curves = [np.loadtxt(path, comments=("#", "@")) for path in paths]
assert len(set(path.resolve() for path in paths)) == 3
assert all(data.shape == (600, 2) and np.isfinite(data).all() for data in curves)
assert all(np.allclose(data[:, 0], curves[0][:, 0]) for data in curves[1:])

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax in axes:
    for i, (data, color, style) in enumerate(
        zip(curves, ("navy", "darkorange", "seagreen"), ("-", "--", ":")), 1
    ):
        ax.plot(data[:, 0], data[:, 1], color=color, linestyle=style,
                label=f"Replicate {i:02d} (2 ns)")
    ax.set_xlabel("Ca-O distance (nm)")
    ax.set_ylabel("RDF g(r)")
    ax.grid(alpha=0.2)
    ax.legend()

axes[0].set_xlim(0, 0.6)
axes[0].set_title("Hydration shells")
axes[1].set_xlim(0.22, 0.29)
axes[1].set_title("First peak: enlarged")
fig.suptitle("Ca IOD - replicates 01, 02 and 03")
fig.tight_layout()
output = folder / "ca_ow_rdf_replicates.png"
assert not output.exists(), output
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"Saved: {output}")
