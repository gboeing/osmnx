#!/bin/bash
set -euo pipefail

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# activate the virtual environment with pre-releases
uv python pin 3.14
uv sync --all-extras --all-groups --upgrade --prerelease=allow
source .venv/bin/activate

# run tests
SKIP=no-commit-to-branch prek run --all-files
pytest --typeguard-packages=osmnx --cov=osmnx --cov-report=term-missing:skip-covered

# restore the environment without pre-releases
uv python pin --rm
uv sync --all-extras --all-groups --upgrade

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
