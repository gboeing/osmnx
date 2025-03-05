#!/bin/bash
set -euo pipefail

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# run the pre-commit hooks for linting/formatting
SKIP=no-commit-to-branch pre-commit run --all-files

# build and validate the package
python -m validate_pyproject ./pyproject.toml
python -m hatch build --clean
python -m twine check --strict ./dist/*

# run the tests and report the test coverage
python -m pytest --verbose --maxfail=1 --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered --numprocesses 3 --dist loadgroup

# build the docs and test that links are alive
python -m sphinx -q -a -E -W --keep-going -b html ./docs/source ./docs/build/html
python -m sphinx -q -a -E -W --keep-going -b linkcheck ./docs/source ./docs/build/linkcheck

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
