#!/bin/bash
coverage run --source osmnx -m pytest --verbose
coverage report -m
