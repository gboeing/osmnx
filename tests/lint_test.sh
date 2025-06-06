#!/bin/bash
set -euo pipefail

# sync the virtual environment
uv sync --all-extras --dev

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# run the pre-commit hooks for linting/formatting
SKIP=no-commit-to-branch uv run pre-commit run --all-files

# build and validate the package
uv run validate-pyproject ./pyproject.toml
uv build
uv run twine check --strict ./dist/*

# run the tests and report the test coverage
uv run pytest --verbose --maxfail=1 --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered --numprocesses=3 --dist=loadgroup

# build the docs and test that links are alive
uv run sphinx-build -q -a -E -W --keep-going -b html ./docs/source ./docs/build/html
uv run sphinx-build -q -a -E -W --keep-going -b linkcheck ./docs/source ./docs/build/linkcheck

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
