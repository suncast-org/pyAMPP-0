#!/usr/bin/env python3

import platform
from setuptools import setup, Extension
from pathlib import Path
import pdb
import sys

common_compile_flags = ['-std=c++11', '-O3', '-fPIC', '-fopenmp']

compile_flags = {
    'Linux': ["-DLINUX", *common_compile_flags],
    'Windows': [*common_compile_flags],
    'Darwin': [*common_compile_flags]
}

link_flags = {
    'Linux': ["-shared", "-fopenmp"],
    'Windows': ["-shared", "-fopenmp"],
    'Darwin': ["-fopenmp"]
}

current_os = platform.system()

if current_os in link_flags:
    extra_link    = link_flags[current_os]
    extra_compile = compile_flags[current_os]
else:
    extra_link    = []
    extra_compile = []

source_dir = Path("./rendergrff/source")
source_files = [str(source_dir / x.name) for x in sorted(source_dir.glob("*.cpp"))]

if current_os != "Windows":
    source_files.remove(str(source_dir / "dllmain.cpp"))

#pdb.set_trace()
source_files.append("pyinit.cpp")

# Define the extension module
render_grff_module = Extension(
    'pyampp0.RenderGRFF',
    sources=source_files,
    extra_compile_args=extra_compile,
    extra_link_args=extra_link,
    include_dirs=[source_dir],
    export_symbols = [],
    language = "c++",
)

setup(
    name="pyampp0",
    ext_modules=[render_grff_module],
)