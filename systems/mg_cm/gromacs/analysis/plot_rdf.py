from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
rdf = np.loadtxt(folder / "mg_ow_rdf.xvg", comments=("#", "@"))
cn = np.loadtxt(folder / "mg_ow_cn.xvg", comments=("#", "@"))

assert rdf.shape[1] == cn.shape[1] == 2
assert np.isfinite(rdf).all() and np.isfinite(cn).all()

cutoff = 0.28
shell = rdf[(rdf[:, 0] >= 0.15) & (rdf[:, 0] <= 0.30)]
peak_r = shell[np.argmax(shell[:, 1]), 0]
coordination = np.interp(cutoff, cn[:, 0], cn[:, 1])

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

axes[0].plot(rdf[:, 0], rdf[:, 1], color="navy")
axes[0].axvline(peak_r, color="green", linestyle=":",
                label=f"First peak: {peak_r:.3f} nm")
axes[0].axhline(1, color="gray", linestyle=":", label="g(r) = 1")
axes[0].set_ylabel("RDF g(r)")

axes[1].plot(cn[:, 0], cn[:, 1], color="darkred")
axes[1].plot(cutoff, coordination, "o", color="black",
             label=f"N({cutoff:.2f} nm) = {coordination:.2f}")
axes[1].set_ylabel("Mean number of water oxygens")
axes[1].set_xlabel("Mg–O distance (nm)")
axes[1].set_ylim(0, np.interp(0.6, cn[:, 0], cn[:, 1]) * 1.05)

for ax in axes:
    ax.axvline(cutoff, color="darkorange", linestyle="--",
               label="Provisional shell cutoff: 0.28 nm")
    ax.set_xlim(0, 0.6)
    ax.grid(alpha=0.2)
    ax.legend()

fig.suptitle("Mg²⁺ CM — hydration structure")
fig.tight_layout()
output = folder / "mg_ow_rdf_cn.png"
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"图片已保存：{output}")
