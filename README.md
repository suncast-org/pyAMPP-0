### What this repository is about?

The aim of this repo is to provide a minimum viable example of [GX Simulator](https://github.com/Gelu-Nita/GX_Simulator) Automated Model Production Pipeline ported line-by-line to Python.

### Installation

Python versions 3.10 and higher are supported. Tested on Windows and Linux. MacOS is not tested and not guaranteed to work yet.

#### Conda env setup

```
conda create -n pyampp_test python=3.10
conda activate pyampp_test
conda install -c conda-forge uv
```

#### Package installation / update

```
pip install uv
uv pip install -U git+https://github.com/vit1-irk/pyAMPP-0
```

### Development installation

```
git clone --recurse-submodules https://github.com/vit1-irk/pyAMPP-0
cd pyAMPP-0
python3 -m pip install -e .
```

### Using with Python

See [test_01.ipynb](./test_01.ipynb).

### Credits

Dependencies:

* <https://github.com/Alexey-Stupishin/pyAMaFiL>
* <https://github.com/kuznetsov-radio/gximagecomputing>
