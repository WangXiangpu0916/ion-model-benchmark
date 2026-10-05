"""Expose pip-installed CUDA DLLs on Windows, then run the unchanged simulator.

Example: python run_gpu.py --campaign pilot --platform CUDA --output results/pilot_cuda
Only this child process's DLL search path is changed; no machine settings change.
"""
import os
from pathlib import Path
import runpy
import sys

if os.name == "nt":
    packages = Path(sys.prefix) / "Lib" / "site-packages"
    directories = [packages / "OpenMM.libs" / "lib", *sorted((packages / "nvidia").glob("*/bin"))]
    os.environ["PATH"] = os.pathsep.join(str(d) for d in directories if d.is_dir()) + os.pathsep + os.environ["PATH"]

runpy.run_path(str(Path(__file__).with_name("simulate.py")), run_name="__main__")
