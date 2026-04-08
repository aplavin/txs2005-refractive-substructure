#!/usr/bin/env bash
set -euo pipefail

REPO="aplavin/txs2005-refractive-substructure"
MSG="Code for: Direct VLBI Detection of Interstellar Turbulence Imprint on a Quasar: TXS 2005+403"

# --- Phase 1: Push clean orphan commit ---
COMMIT=$(git commit-tree HEAD^{tree} -m "$MSG")
git push origin "$COMMIT":refs/heads/master --force
echo "Published as $COMMIT"

# --- Opt-in: Create GitHub release for Zenodo ---
if [[ "${1:-}" == "--release" ]]; then
    VERSION="${2:?Usage: $0 --release <version>}"

    git tag "$VERSION" "$COMMIT"
    git push origin "$VERSION"

    gh release create "$VERSION" \
        --repo "$REPO" \
        --target "$COMMIT" \
        --title "$VERSION: $MSG"

    echo "Release: https://github.com/$REPO/releases/tag/$VERSION"
fi
