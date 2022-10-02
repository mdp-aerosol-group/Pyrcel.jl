using Conda

Conda.pip_interop(true)
Conda.add(["assimulo", "numpy", "cython", "numba", "future", "pandas"])
Conda.pip("install", "pyrcel")
