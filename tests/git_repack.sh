#!/bin/bash
docker image prune -af && docker system prune -af && docker volume prune -af
git remote prune origin
git repack -a -d -f --depth=250 --window=250
git prune-packed
git reflog expire --expire=1.month.ago
git gc --aggressive
