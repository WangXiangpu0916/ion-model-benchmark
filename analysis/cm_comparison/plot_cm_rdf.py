from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
project = folder.parents[1]
fig, ax = plt.subplots(figsize=(9, 5))

for ion, label, color in [
    ("mg", "Mg²⁺ CM", "navy"),
    ("ca", "Ca²⁺ CM", "darkorange"),
]:
    source = (
        project / "systems" / f"{ion}_cm"
        / "gromacs" / "analysis" / f"{ion}_ow_rdf.xvg"
    )
    data = np.loadtxt(source, comments=("#", "@"))
    assert data.ndim == 2 and data.shape[1] == 2
    assert np.isfinite(data).all()
    assert np.all(np.diff(data[:, 0]) > 0)
    ax.plot(data[:, 0], data[:, 1], label=label, color=color)

ax.axhline(1, color="gray", linestyle=":", linewidth=1)
ax.set_xlim(0.15, 0.60)
ax.set_ylim(bottom=0)
ax.set_xlabel("Ion–water oxygen distance (nm)")
ax.set_ylabel("RDF g(r)")
ax.set_title("CM pilots: Mg²⁺ vs Ca²⁺ (2 ns each)")
ax.grid(alpha=0.2)
ax.legend()
fig.tight_layout()

output = folder / "cm_rdf_comparison.png"
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"Saved: {output}")
