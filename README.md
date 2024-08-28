### What this repository is about?

The aim of this repo is to provide a minimum viable example of [GX Simulator](https://github.com/Gelu-Nita/GX_Simulator) Automated Model Production Pipeline ported line-by-line to Python.

```bash
# without radio source modeling support
python3 -m pip install git+https://github.com/vit1-irk/pyAMPP-0

# with radio source modeling support

git clone --recurse-submodules https://github.com/vit1-irk/pyAMPP-0
cd pyAMPP-0
python3 -m pip install .
```

### Using with Python

See [test_01.ipynb](./test_01.ipynb).

### Credits

Upstream submodules:

* <https://github.com/Alexey-Stupishin/pyAMaFiL>
* <https://github.com/kuznetsov-radio/gximagecomputing>
