#!/usr/bin/env python
"""Generate requirements-all.txt for creating a conda environment."""

from itertools import chain
from pathlib import Path

from tomllib import load

path_pyproject = "../pyproject.toml"
path_extras = "./requirements-extras.txt"
path_output = "./requirements-all.txt"

with Path(path_pyproject).open("rb") as f:
    pyproject = load(f)
with Path(path_extras).open("rt") as f:
    extras = f.read().splitlines()

deps = pyproject["project"]["dependencies"]
devs = pyproject["dependency-groups"]["dev"]
opts = list(chain.from_iterable(pyproject["project"]["optional-dependencies"].values()))
reqs = sorted(set(deps + opts + devs + extras))

with Path(path_output).open("w") as f:
    f.writelines("\n".join(reqs))
