#!/bin/bash

# exit on error
set -e

# delete temp files and folders
rm -r -f .coverage .pytest_cache .temp ./docs/_build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete

# check if imports are organized properly
isort . --check-only

# check if code is formatted properly
black . --check --diff

# lint the code
flake8 .

# lint the docstrings
pydocstyle .

# build the docs
# make -C ./docs html
# python -m sphinx -b linkcheck docs/ docs/_build/linkcheck

# run the tests
coverage run --source ./osmnx --module pytest --verbose

# report the test coverage
coverage report -m

# delete temp files and folders
rm -r -f .coverage .pytest_cache .temp ./docs/_build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
