#!/bin/bash
set -euo pipefail

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__

# activate the virtual environment with minimal versions
uv python pin 3.11
uv sync --all-extras --group test --upgrade --resolution=lowest-direct
source .venv/bin/activate

# run tests
python ./tests/verify_min_deps.py
pytest -W ignore

# restore the environment with current release versions
uv python pin --rm
uv sync --all-extras --all-groups --upgrade

# delete temp files and folders
rm -r -f ./.coverage* ./.pytest_cache ./.temp ./dist ./docs/build ./*/__pycache__
