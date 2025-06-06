#!/usr/bin/env python
# /// script
# dependencies = [
#     "packaging",
#     "tomli",
# ]
# ///
"""Generate requirements-min.txt for testing minimum dependency versions."""

from itertools import chain
from pathlib import Path

from packaging.requirements import Requirement
from tomli import load as toml_load

# load the pyproject.toml file
with Path("./pyproject.toml").open("rb") as f:
    pyproject = toml_load(f)

# extract/pin dependencies + optionals, and extract dev group from pyproject
deps = [Requirement(d) for d in pyproject["project"]["dependencies"]]
optionals = pyproject["project"]["optional-dependencies"].values()
deps.extend({Requirement(o) for o in chain.from_iterable(optionals)})
reqs = [f"{dep.name}=={next(iter(dep.specifier)).version}.*" for dep in deps]
reqs.extend(pyproject["dependency-groups"]["dev"])

# write the requirements as to file on disk
with Path("./requirements-min.txt").open("w") as f:
    f.writelines("\n".join(sorted(reqs)))
