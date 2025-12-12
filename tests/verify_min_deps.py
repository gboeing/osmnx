#!/usr/bin/env python
"""Verify that installed dependencies match minimum dependency versions."""

from importlib.metadata import version as metadata_version
from itertools import chain
from pathlib import Path
from tomllib import load as toml_load

from packaging.requirements import Requirement

# load the pyproject.toml file
with Path("./pyproject.toml").open("rb") as f:
    pyproject = toml_load(f)

# extract and pin all required + optional dependencies from pyproject
deps = [Requirement(d) for d in pyproject["project"]["dependencies"]]
opts = [v for k, v in pyproject["project"]["optional-dependencies"].items() if k != "all"]
deps.extend({Requirement(o) for o in chain.from_iterable(opts)})
requirements = {dep.name: next(iter(dep.specifier)).version for dep in deps}
requirements = dict(sorted(requirements.items()))

# check that installed versions match minimum versions
ok_msg = ""
err_msg = ""
for package, required_version in requirements.items():
    installed_version = metadata_version(package)
    if installed_version.startswith(required_version):
        ok_msg += f"\nExpected {package} {required_version}, matches {installed_version}."
    else:
        err_msg += f"\nExpected {package} {required_version}, found {installed_version}."

# print ok message or raise error with error message
if err_msg == "":
    ok_msg = "Installed dependencies match minimum dependency versions." + ok_msg
    print(ok_msg)  # noqa: T201
else:
    err_msg = "Installed dependencies do not match minimum dependency versions." + err_msg
    raise ImportError(err_msg)
