#!/bin/bash
pylint -E ./osmnx/*
coverage run --source osmnx --module pytest --verbose
coverage report -m
