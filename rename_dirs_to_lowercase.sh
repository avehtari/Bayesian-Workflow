#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

dirs=(
  "Bioassay"
  "Coronavirus"
  "Declining_exponentials"
  "Dogs"
  "Golf"
  "Movies"
  "Multiple_choice"
  "Nabiximols"
  "Park_rule"
  "Planetary_motion"
  "Roaches"
  "Sleepstudy"
  "Timeseries"
  "Worldcup"
)

for d in "${dirs[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "Skipping: '$d' (not found in repo root)"
    continue
  fi

  lower="$(printf '%s' "$d" | tr '[:upper:]' '[:lower:]')"

  if [[ "$d" == "$lower" ]]; then
    echo "Skipping: '$d' (already lowercase)"
    continue
  fi

  tmp="__tmp__${lower}__$$"
  echo "Renaming: '$d' -> '$lower'"
  git mv -f "$d" "$tmp"
  git mv -f "$tmp" "$lower"
done

echo
echo "Done. Review with: git status"
echo "Then commit: git commit -m \"Rename directories to lowercase\""
