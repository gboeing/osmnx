#!/bin/bash
ruff check . --fix-only --no-cache
black .
