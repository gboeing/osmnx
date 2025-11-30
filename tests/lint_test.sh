#!/bin/bash
set -euo pipefail

# activate the virtual environment
source .venv/bin/activate

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# run the pre-commit hooks for linting/formatting
SKIP=no-commit-to-branch prek run --all-files

# build and validate the package
uv build
twine check --strict ./dist/*
validate-pyproject ./pyproject.toml

# run the tests and report the test coverage
pytest --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered

# build the docs and test that links are alive
sphinx-build -q -a -E -W --keep-going -b html ./docs/source ./docs/build/html
sphinx-build -q -a -E -W --keep-going -b linkcheck ./docs/source ./docs/build/linkcheck

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
