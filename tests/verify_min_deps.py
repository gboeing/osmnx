#!/usr/bin/env python
"""Verify installed dependencies match minimum dependency versions."""

from importlib.metadata import version as metadata_version
from itertools import chain
from pathlib import Path

from packaging.requirements import Requirement
from tomli import load as toml_load

# load the pyproject.toml file
with Path("./pyproject.toml").open("rb") as f:
    pyproject = toml_load(f)

# extract/pin dependencies + optionals + dev group from pyproject
deps = [Requirement(d) for d in pyproject["project"]["dependencies"]]
optionals = pyproject["project"]["optional-dependencies"].values()
deps.extend({Requirement(o) for o in chain.from_iterable(optionals)})
deps.extend(Requirement(d) for d in pyproject["dependency-groups"]["dev"])
requirements = {dep.name: next(iter(dep.specifier)).version for dep in deps}
requirements = dict(sorted(requirements.items()))

message = ""
for package, required_version in requirements.items():
    installed_version = metadata_version(package)
    if not installed_version.startswith(required_version):
        message += f"Expected {package} {required_version}, found {installed_version}. "

if message != "":
    raise ImportError(message)
