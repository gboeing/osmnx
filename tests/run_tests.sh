#!/bin/bash
coverage run --source osmnx --module pytest --verbose
coverage report -m
