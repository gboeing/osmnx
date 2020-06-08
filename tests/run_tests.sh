#!/bin/bash
set -e
flake8
pydocstyle
coverage run --source osmnx --module pytest --verbose
coverage report -m
