#!/bin/bash
black --line-length 100 ./setup.py
black --line-length 100 ./docs/source/conf.py
black --line-length 100 ./tests/test_osmnx.py
black --line-length 100 ./osmnx/*
