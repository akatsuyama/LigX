#実行は　python setup.py sdistと python setup.py bdist_wheel

from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
import shutil, re

exec(open('LigX/_version.py').read())

NAME = "ligx"
DESCRIPTION = "LigX: ligand generation in protein binding site."
AUTHOR = "Akira Katsuyama"
AUTHOR_EMAIL = "katsuyama@pharm.hokudai.ac.jp"
URL = "https://github.com/akatsuyama/LigX"
DOWNLOAD_URL = "https://github.com/akatsuyama/LigX"
LICENSE = "MIT"
PYTHON_REQUIRES = ">=3.9"


VERSION = __version__

def read_long_description():
    readme_md = Path("README.md")
    if readme_md.exists():
        return readme_md.read_text(encoding="utf-8"), "text/markdown"
    return (DESCRIPTION, "text/plain")

long_description, long_description_content_type = read_long_description()


ASSET_DIRS = ["fragment", "rollout_fragment"]

class build_py(_build_py):
    def run(self):
        super().run()
        # copy to build/lib/ligx/<dir>/
        for src_name in ASSET_DIRS:
            src = Path(src_name)
            if not src.exists():
                continue
            target = Path(self.build_lib) / "LigX" / src_name
            target.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, target, dirs_exist_ok=True)  # 3.8+



INSTALL_REQUIRES = [
    "numpy>=1.20",
]


PACKAGES = [
    'LigX'
]

CLASSIFIERS = [
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Operating System :: OS Independent",
]


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=Path("README.md").read_text(encoding="utf-8") if Path("README.md").exists() else DESCRIPTION,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    license=LICENSE,
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    packages=PACKAGES,
    include_package_data=True,  # MANIFEST.in
    package_data={
        "LigX": [
            "fragment/*", "fragment/**/*",
            "rollout_fragment/*", "rollout_fragment/**/*",
            "clean.sh"
        ],
    },
    cmdclass={"build_py": build_py},
    classifiers=CLASSIFIERS,
    keywords=["compound-generation", "drug-discovery", "structure-based drug design", "medicinal-chemistry"],
)
