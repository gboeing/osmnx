#!/bin/bash

# exit on error
set -e

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist/ ./docs/build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
find . -type f -name ".coverage*" -delete

# lint
#ruff check . --preview
pre-commit run --all-files

# test building and validating the package
hatch build --clean
twine check --strict ./dist/*

# build the docs
make -C ./docs html SPHINXOPTS="-W --keep-going"
#python -m sphinx -b linkcheck ./docs/source ./docs/build/linkcheck

# run the tests and report the test coverage
pytest --maxfail=1 --typeguard-packages=osmnx --cov=./osmnx --cov-report=term-missing --verbose

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist/ ./docs/build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
find . -type f -name ".coverage*" -delete
