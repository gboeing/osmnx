#!/bin/bash
set -e

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist ./docs/build ".coverage*" "./*/__pycache__"

# run the pre-commit hooks for linting/formatting
pre-commit run --all-files

# build and validate the package
python -m validate_pyproject ./pyproject.toml
python -m hatch build --clean
python -m twine check --strict ./dist/*

# build the docs and test that links are alive
python -m sphinx -E -W --keep-going -b html ./docs/source ./docs/build/html
python -m sphinx -E -W --keep-going -b linkcheck ./docs/source ./docs/build/linkcheck

# run the tests and report the test coverage
python -m pytest --verbose --maxfail=1 --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist ./docs/build ".coverage*" "./*/__pycache__"
