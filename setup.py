
'''
patch: renormhnet
'''
import os.path
import codecs
from setuptools import setup, find_packages


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="renormhnet",
    version=get_version("renormhnet/__init__.py"),
    description="RenormHNet: Renormalization of Hamiltonian network currently supports Heisenberg exchange (Jx) Hamiltonians.",
    download_url="https://github.com/DHKiem/RenormHNet",
    author="He Li (DeepH-pack); Do Hoon Kiem (RenormHNet)",
    python_requires=">=3.9",
    packages=find_packages(),
    package_dir={'renormhnet': 'renormhnet'},
    package_data={'': ['*.jl', '*.ini', 'periodic_table.json']},
    entry_points={
        "console_scripts": [
            "deeph-preprocess = renormhnet.scripts.preprocess:main",
            "deeph-train = renormhnet.scripts.train:main",
            "deeph-evaluate = renormhnet.scripts.evaluate:main",
            "deeph-inference = renormhnet.scripts.inference:main",
            "inference-jx = renormhnet.scripts.inference_jx:main",
            "jxbinder = renormhnet.scripts.jxbinder:main",
        ]
    },
    install_requires=[
        "numpy",
        "scipy",
        "torch>=1.9",
        "torch_geometric>=1.7.2",
        "e3nn>=0.3.5, <=0.4.4",
        "h5py",
        "pymatgen",
        "pathos",
        "psutil",
        "tqdm",
        "tensorboard",
        "pandas",
        "matplotlib",
        "tomli; python_version < '3.11'",
    ],
    license="LGPL-3.0",
    license_files="LICENSE",
    zip_safe=False,
)
