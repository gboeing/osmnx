#!/bin/bash
set -e  # exit on error

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist/ ./docs/build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
find . -type f -name ".coverage*" -delete

# lint
pre-commit run --all-files

# build the docs
make -C ./docs html SPHINXOPTS="-E -W --keep-going"

# run the tests and report the test coverage
pytest --verbose --maxfail=1 --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered

# delete temp files and folders
rm -r -f .pytest_cache .temp ./dist/ ./docs/build osmnx/__pycache__ tests/__pycache__
find . -type f -name "*.vrt" -delete
find . -type f -name ".coverage*" -delete
