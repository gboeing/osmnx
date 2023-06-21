#!/bin/bash
ruff check . --fix-only
black .
