#!/bin/bash
isort . -rc -y
black . --line-length 100
