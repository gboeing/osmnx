#!/bin/bash
set -euo pipefail

CACHE_DIR="tests/.cache"
CACHE_READY_FILE="${CACHE_DIR}/.cache_ready"
USE_OSMNX_CACHE="${USE_OSMNX_CACHE:-true}"

# use four workers only with a ready cache and cache enabled
WORKERS=1
if [[ "$USE_OSMNX_CACHE" == "true" && -f "$CACHE_READY_FILE" ]]; then
    WORKERS=4
fi

echo "Running tests with ${WORKERS} worker(s); USE_OSMNX_CACHE=${USE_OSMNX_CACHE}."
pytest --numprocesses="$WORKERS" --dist=loadgroup "$@"

# mark the cache ready after every successful cache-enabled run
if [[ "$USE_OSMNX_CACHE" == "true" ]]; then
    touch "$CACHE_READY_FILE"
fi
