#!/bin/bash
set -euo pipefail

# clean the caches of pip and uv
pip cache purge
uv cache clean

# deactivate conda and clean its cache
eval "$(conda shell.bash hook)"
conda deactivate
conda clean --all --yes

# remove unused docker images, data, and local volumes
docker image prune -af && docker system prune -af && docker volume prune -af

# prune, repack, and garbage collect git
git remote prune origin
git repack -a -d -f --depth=250 --window=250
git prune-packed
git reflog expire --expire=1.month.ago
git gc --aggressive
