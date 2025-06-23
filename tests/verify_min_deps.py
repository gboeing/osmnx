#!/usr/bin/env python
"""Verify that installed dependencies match minimum dependency versions."""

from importlib.metadata import version as metadata_version
from itertools import chain
from pathlib import Path

from packaging.requirements import Requirement
from tomli import load as toml_load

# load the pyproject.toml file
with Path("./pyproject.toml").open("rb") as f:
    pyproject = toml_load(f)

# extract/pin dependencies + all extras + all group from pyproject
deps = [Requirement(d) for d in pyproject["project"]["dependencies"]]
opts = pyproject["project"]["optional-dependencies"].values()
deps.extend({Requirement(o) for o in chain.from_iterable(opts) if not o.startswith("osmnx")})
deps.extend(Requirement(d) for g in pyproject["dependency-groups"].values() for d in g)

requirements = {dep.name: next(iter(dep.specifier)).version for dep in deps}
requirements = dict(sorted(requirements.items()))

ok_msg = "Installed dependencies match minimum dependency versions."
err_msg = ""
for package, required_version in requirements.items():
    installed_version = metadata_version(package)
    if installed_version.startswith(required_version):
        ok_msg += f"\nExpected {package} {required_version}, matches {installed_version}."
    else:
        err_msg += f"\nExpected {package} {required_version}, found {installed_version}."

if err_msg != "":
    err_msg = "Installed dependencies do not match minimum dependency versions." + err_msg
    raise ImportError(err_msg)

print(ok_msg)  # noqa: T201
