#!/bin/bash
coverage run --source osmnx --module pytest --verbose
coverage report -m
flake8 .
pydocstyle .
