#!/bin/bash
set -e
flake8
pydocstyle
coverage run --source osmnx --module pytest --verbose
coverage report -m
rm -r -f .coverage .pytest_cache .temp
