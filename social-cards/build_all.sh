#!/usr/bin/env bash
# ── Build every case-study social card from cards.tsv ──
# Usage: social-cards/build_all.sh [docs_dir]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DOCS_DIR="${1:-}"

grep -vE '^\s*(#|$)' "$HERE/cards.tsv" | while IFS=$'\t' read -r page fig; do
  [ -n "$page" ] || continue
  "$HERE/make_card.sh" "$page" "$fig" ${DOCS_DIR:+"$DOCS_DIR"} || echo "  ^ FAILED: $page" >&2
done
