### What this repository is about?

The aim of this repo is to provide a minimum viable example of [GX Simulator](https://github.com/Gelu-Nita/GX_Simulator) Automated Model Production Pipeline ported line-by-line to Python.

Clone with all of the submodules 

```bash
git clone https://github.com/suncast-org/pyAMPP-0
cd pyAMPP-0
git submodule update --init --remote --recursive
```

and run `example.ipynb`

Don't forget to compile Linux version of library in `nlfff` directory!

### Does it work?

It's not yet fully working, you still need IDL-produced .SAV files for potential field and photosphere bounded boxes.

### Credits

Upstream submodules:

* <https://github.com/Alexey-Stupishin/Magnetic-Field_Library>
* <https://github.com/kuznetsov-radio/gximagecomputing>