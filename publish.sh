#!/usr/bin/env bash
set -euo pipefail

MSG="Code for: Direct VLBI Detection of Interstellar Turbulence Imprint on a Quasar: TXS 2005+403"
COMMIT=$(git commit-tree HEAD^{tree} -m "$MSG")
git push origin "$COMMIT":refs/heads/master --force
echo "Published as $COMMIT"
