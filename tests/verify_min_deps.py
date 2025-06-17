#!/usr/bin/env python
# /// script
# dependencies = [
#     "packaging",
#     "tomli",
# ]
# ///
"""Verify installed dependencies match minimum dependency versions."""

from importlib.metadata import version
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
reqs = {dep.name: next(iter(dep.specifier)).version for dep in deps}
reqs = dict(sorted(reqs.items()))

wrong_versions = []
for pkg, v in reqs.items():
    installed_version = version(pkg)
    if not installed_version.startswith(v):
        wrong_versions.append((pkg, v, installed_version))

if len(wrong_versions) > 0:
    msg = ""
    for pkg, v, installed_version in wrong_versions:
        msg += f"Expected {pkg} {v}, found {installed_version}. "
    raise ImportError(msg)
