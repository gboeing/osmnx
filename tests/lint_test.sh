#!/bin/bash
set -euo pipefail

# activate the virtual environment
source .venv/bin/activate

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# run the pre-commit hooks for linting/formatting
SKIP=no-commit-to-branch pre-commit run --all-files

# build and validate the package
uv build
uvx twine check --strict ./dist/*
uvx --from=validate-pyproject[all] validate-pyproject ./pyproject.toml

# run the tests and report the test coverage
pytest --verbose --maxfail=1 --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered --numprocesses=3 --dist=loadgroup

# build the docs and test that links are alive
sphinx-build -q -a -E -W --keep-going -b html ./docs/source ./docs/build/html
sphinx-build -q -a -E -W --keep-going -b linkcheck ./docs/source ./docs/build/linkcheck

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
