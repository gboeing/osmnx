#!/bin/bash

# exit on error
set -e

# delete temp files and folders
rm -r -f .coverage .pytest_cache .temp ./dist/ ./docs/_build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete

# test building and validating the package
hatch build --clean
twine check --strict ./dist/*

# lint
pre-commit run --all-files

# build the docs
# make -C ./docs html
# python -m sphinx -b linkcheck docs/ docs/_build/linkcheck

# run the tests and report the test coverage
pytest --cov=./osmnx --cov-report=term-missing --verbose

# delete temp files and folders
rm -r -f .coverage .pytest_cache .temp ./dist/ ./docs/_build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
