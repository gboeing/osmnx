#!/bin/bash
set -euo pipefail

CACHE_DIR="tests/.cache"
CACHE_READY_FILE="${CACHE_DIR}/.cache_ready"
USE_PERSISTENT_CACHE="${USE_PERSISTENT_CACHE:-true}"

# use four workers only with a persistent cache that's been marked ready
WORKERS=1
if [[ "$USE_PERSISTENT_CACHE" == "true" && -f "$CACHE_READY_FILE" ]]; then
    WORKERS=4
fi

# if we're not using the persistent cache, verify the temp cache is cleared
if [[ "$USE_PERSISTENT_CACHE" != "true" ]]; then
    rm -r -f ./tests/.temp/cache/
fi

echo "Running tests with ${WORKERS} worker(s); USE_PERSISTENT_CACHE=${USE_PERSISTENT_CACHE}."
pytest --numprocesses="$WORKERS" --dist=loadgroup "$@"

# mark the persistent cache ready after every successful persistent-cache run
if [[ "$USE_PERSISTENT_CACHE" == "true" ]]; then
    touch "$CACHE_READY_FILE"
fi
