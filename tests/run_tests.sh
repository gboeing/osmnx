#!/bin/bash
flake8
coverage run --source osmnx --module pytest --verbose
coverage report -m
