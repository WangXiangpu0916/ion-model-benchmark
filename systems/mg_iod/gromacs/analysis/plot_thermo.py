from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
data = np.loadtxt(folder / "md_thermo.xvg", comments=("#", "@"))

assert data.shape[1] == 6, "应该有时间列和五列数据"
assert np.isfinite(data).all(), "数据中存在无效数值"
assert np.allclose(np.diff(data[:, 0]), 0.1), "请检查数据的时间间隔"

time_ns = data[:, 0] / 1000

panels = [
    (data[:, 1], "Potential energy (kJ/mol)", None),
    (data[:, 2], "Temperature (K)", 298),
    (data[:, 3], "Pressure (bar)", 1),
    (data[:, 4], "Volume (nm³)", None),
    (data[:, 5] / 1000, "Density (g/mL)", None),
]

fig, axes = plt.subplots(5, 1, figsize=(10, 12), sharex=True)

for ax, (values, label, target) in zip(axes, panels):
    smoothed = pd.Series(values).rolling(501, center=True).mean()
    ax.plot(time_ns, values, color="steelblue", alpha=0.3,
            linewidth=0.5, label="Recorded values")
    ax.plot(time_ns, smoothed, color="navy",
            linewidth=1.5, label="About 50 ps moving average")
    if target is not None:
        ax.axhline(target, color="darkorange", linestyle="--",
                   linewidth=1, label="Target")
    ax.set_ylabel(label)
    ax.grid(alpha=0.2)

axes[0].legend(loc="best")
axes[1].legend(loc="best")
axes[-1].set_xlabel("Time (ns)")
fig.suptitle("Mg²⁺ IOD — replicate 01 - 2 ns sampling")
fig.tight_layout()

output = folder / "md_thermo.png"
assert not output.exists(), output
fig.savefig(output, dpi=180)
plt.close(fig)
print(f"图片已保存：{output}")
